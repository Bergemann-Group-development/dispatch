import re
from pathlib import Path
from numbers import Integral, Real

# --- regexes ---
_QUOTE_RE = re.compile(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"")
_GROUP_RE = re.compile(r"&(\w+)(.*?)\/", flags=re.S)
_KEY_RE = re.compile(r"(\w+)\s*=")
_REPEAT_RE = re.compile(r"^(\d+)\*(.+)$")

def _unquote(orig, strip=True):
    """Remove quotes and unescape doubled quotes."""
    if orig.startswith("'") and orig.endswith("'"):
        s = orig[1:-1].replace("''", "'")
    elif orig.startswith('"') and orig.endswith('"'):
        s = orig[1:-1].replace('""', '"')
    else:
        s = orig
    return s.strip() if strip else s

def _parse_value(tok):
    low = tok.lower()
    if low in ('.true.', 't', 'true'):
        return True
    if low in ('.false.', 'f', 'false'):
        return False
    try:
        return int(tok)
    except ValueError:
        try:
            return float(tok.replace('d', 'e'))  # Fortran 1.0d0 → 1.0e0
        except ValueError:
            return tok

def _format_value(val):
    """Format Python value back to Fortran namelist syntax."""
    if isinstance(val, bool):
        return ".true." if val else ".false."
    elif isinstance(val, (Integral, Real)):
        return str(val)
    elif isinstance(val, str):
        # Always single-quote strings; escape by doubling
        return f"'{val}'"
    elif isinstance(val, Path):
        return f"'{str(val)}'"
    elif val is None:
        return ''
    else:
        raise TypeError(f"Unsupported type for namelist output: {type(val)}")


class NameList(dict):
    """Simple Fortran Namelist parser/wrapper."""

    def __init__(self, data=None):
        super().__init__(data or {})

    # --- Constructors ---
    @classmethod
    def from_string(cls, text, strip_strings=True):
        return cls(_parse_namelist(text, strip_strings))

    @classmethod
    def from_file(cls, filename, strip_strings=True, encoding="utf-8"):
        text = Path(filename).read_text(encoding=encoding)
        return cls.from_string(text, strip_strings=strip_strings)

    # --- Convenience ---
    def groups(self):
        """Return list of group names."""
        return list(self.keys())

    def __getitem__(self, key):
        return super().__getitem__(key.lower())

    def __setitem__(self, key, value):
        return super().__setitem__(key.lower(), value)

    def __contains__(self, key):
        return super().__contains__(key.lower())

    # --- Output ---
    def to_string(self):
        """Return Fortran namelist string."""
        lines = []
        for gname, gdict in self.items():
            lines.append(f"&{gname}")
            for key, val in gdict.items():
                if isinstance(val, list):
                    vals = " ".join(_format_value(v) for v in val)
                else:
                    vals = _format_value(val)
                lines.append(f"  {key} = {vals}")
            lines.append("/")
        return "\n".join(lines)

    def to_file(self, filename, encoding="utf-8"):
        """Write namelist to file."""
        Path(filename).write_text(self.to_string(), encoding=encoding)

    def get_item(self, namelist, key, fallback=None):
        """Get item from group with fallback."""
        return self.get(namelist, {}).get(key, fallback)



# --- internal parser ---
def _parse_namelist(text, strip_strings=True):
    # 1) Extract quotes -> placeholders
    quotes = []
    def _quote_replacer(m):
        idx = len(quotes)
        quote = m.group(0)
        quote = quote.replace('\n', '')  # remove newlines inside quotes
        quotes.append(quote)
        return f"__STR{idx}__"

    text_with_ph = _QUOTE_RE.sub(_quote_replacer, text)

    # 2) Remove comments
    text_nocom = re.sub(r'!.*$', '', text_with_ph, flags=re.M)

    # 3) Normalize
    text_norm = text_nocom.replace(',', ' ')
    placeholder_map = {f"__STR{i}__": quotes[i] for i in range(len(quotes))}

    result = {}

    # 4) Extract &name ... /
    for gname, body in _GROUP_RE.findall(text_norm):
        gname = gname.lower()
        gdict = {}

        # 5) Find all key=value
        key_matches = list(_KEY_RE.finditer(body))
        for idx, m in enumerate(key_matches):
            key = m.group(1).lower()
            start = m.end()
            end = key_matches[idx+1].start() if idx+1 < len(key_matches) else len(body)
            raw_val = body[start:end].strip()
            tokens = raw_val.split()

            parsed_vals = []
            for tok in tokens:
                if tok in placeholder_map:
                    parsed_vals.append(_unquote(placeholder_map[tok], strip=strip_strings))
                    continue

                mrep = _REPEAT_RE.match(tok)
                if mrep:
                    n = int(mrep.group(1))
                    val = _parse_value(mrep.group(2))
                    if val in placeholder_map:
                        val = _unquote(placeholder_map[val], strip=strip_strings)
                    parsed_vals.extend([val] * n if not isinstance(val, list) else val * n)
                else:
                    parsed_vals.append(_parse_value(tok))

            # simplify singleton
            value = parsed_vals if len(parsed_vals) > 1 else (parsed_vals[0] if parsed_vals else None)

            # merge repeated keys
            if key in gdict:
                if isinstance(gdict[key], list):
                    if isinstance(value, list):
                        gdict[key].extend(value)
                    else:
                        gdict[key].append(value)
                else:
                    gdict[key] = [gdict[key]] + (value if isinstance(value, list) else [value])
            else:
                gdict[key] = value

        result[gname] = gdict

    return result