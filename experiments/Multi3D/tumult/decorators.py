import functools
import numpy as np
import inspect
import operator

def broadcast_arguments(obj='self', numeric=True, **explicit_dims):
    def decorator(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            # Extract actual call arguments (positional override + kwargs)
            arguments = {k: v for i, (k, v) in enumerate(bound.arguments.items()) if i < len(args)}
            arguments.update(kwargs)

            # Handle possible method call
            self = arguments.pop(obj, None)
            is_obj = self is not None

            # ---------- STEP 1: preprocess input arguments ----------
            info = {}
            for name, val in arguments.items():
                arr = np.asarray(val)
                defdim = explicit_dims.get(name, 0)
                extradim = max(arr.ndim - defdim, 0)
                tail  = arr.shape[extradim:]          # dims preserved for func

                info[name] = (arr, extradim, tail)

            # ---------- STEP 2: compute broadcast shape for extra dims ----------
            extra_shapes = [arr.shape[:extradim] for arr, extradim, tail in info.values()]
            target_extra = np.broadcast_shapes(*extra_shapes) if extra_shapes else ()

            # ---------- STEP 3: broadcast *every* argument (even scalars) ----------
            bcast = {}
            for name, (arr, extradim, tail) in info.items():
                # Everything becomes shape = target_extra + tail
                full_shape = target_extra + tail
                bcast[name] = np.broadcast_to(arr, full_shape)

            # ---------- STEP 4: iterate over broadcasted grid ----------
            results = np.empty(target_extra, dtype=object)

            for idx in np.ndindex(target_extra):
                new_kwargs = {name: arr[idx] for name, arr in bcast.items()}
                results[idx] = func(self, **new_kwargs) if is_obj else func(**new_kwargs)

            # ---------- STEP 5: numeric conversion helper ----------
            def convert_to_numeric(arr):
                if not numeric:
                    return arr
                try:
                    return np.array(arr.tolist())
                except Exception:
                    return arr

            # ---------- STEP 6: handle result shape & tuple-return ----------
            example = results.flat[0]

            # case: function returns tuple → split into per-component arrays
            if isinstance(example, tuple):
                ntup = len(example)
                outs = [np.empty(target_extra, dtype=object) for _ in range(ntup)]

                for idx in np.ndindex(target_extra):
                    tpl = results[idx]
                    for i in range(ntup):
                        outs[i][idx] = tpl[i]

                outs = [convert_to_numeric(arr) for arr in outs]
                return outs

            # case: single output
            final = convert_to_numeric(results)

            # reduce scalar arrays
            if final.size == 1:
                return final.item()

            return final
        return wrapper
    return decorator


def auto_ops(attr):
    def decorator(cls):
        opmap = {
            "__add__": operator.add,
            "__sub__": operator.sub,
            "__mul__": operator.mul,
            "__truediv__": operator.truediv,
            "__radd__": operator.add,
            "__rsub__": lambda a, b: operator.sub(b, a),
            "__rmul__": operator.mul,
            "__rtruediv__": lambda a, b: operator.truediv(b, a),
        }
        for name, op in opmap.items():
            def make_method(op):
                def method(self, other):
                    a = getattr(self, attr)
                    b = getattr(other, attr) if isinstance(other, cls) else other
                    return cls(op(a, b), base=self)
                return method
            setattr(cls, name, make_method(op))
        return cls
    return decorator