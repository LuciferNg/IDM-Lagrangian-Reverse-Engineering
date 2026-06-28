# CLASS with IDM Fluid Extension

CLASS (Cosmic Linear Anisotropy Solving System) with the IDM sync dip
fluid equation of state patched into `source/background.c`.

## IDM Modification

The IDM fluid implements:
```
f(z) = 1 - eps*(z/zc)*exp(-z/zc)
w(z) = -1 + (1+z)/(3*zc) * f'(z)/f(z)
```

Parameters in `.ini`:
- `fluid_equation_of_state=IDM`
- `eps_idm=0.1545`
- `zc_idm=0.6`

## Build

```bash
cd class
make clean
make class
```

Requires: gcc, libgsl (GNU Scientific Library).

## Source Files

| Directory | Content |
|:----------|:--------|
| `source/*.c` | CLASS core modules with IDM patches |
| `include/*.h` | Header files |
| `external/` | External dependency headers |
| `tools/` | C utility libraries |
| `main/` | Entry point |
| `python/` | Python wrapper (classy) |

## Test Configs

- `test_lcdm.ini` — Standard ΛCDM (Planck-like parameters)
- `test_idm.ini` — IDM with ε=0.1545, zc=0.6