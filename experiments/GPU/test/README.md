# GCC

## General
- Offloaded code does not work with values from inside user defined types
    - Values and pointers must be placed in seperate variables

## 9.2.1
- Kernels can not be in multiple files and linked together

## 9.3.0
- Kernels can now be in multiple files and be linked together
- Pointers can not be created in offloaded functions as it will give an out-of-bounds memory error
- teams directives must be placed close to the target directives if used together
