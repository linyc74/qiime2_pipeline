The following LEfSe error is caused by newer versions of `rpy2`.

```
AttributeError: 'NoneType' object has no attribute 'rownames'
```

`rpy2` needs to be downgraded to `3.5.10`.

```bash
pip install rpy2==3.5.10
```

Also another LEfSe error caused by `rpy2`. The solution is the same, install version `3.5.10`.

```
Error: package or namespace load failed for ‘methods’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘datasets’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘utils’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘grDevices’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘graphics’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘stats’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Error: package or namespace load failed for ‘methods’:
 .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
During startup - Warning messages:
1: package "methods" in options("defaultPackages") was not found
2: package ‘datasets’ in options("defaultPackages") was not found
3: package ‘utils’ in options("defaultPackages") was not found
4: package ‘grDevices’ in options("defaultPackages") was not found
5: package ‘graphics’ in options("defaultPackages") was not found
6: package ‘stats’ in options("defaultPackages") was not found
7: package ‘methods’ in options("defaultPackages") was not found
Error: .onLoad failed in loadNamespace() for 'methods', details:
  call: assign(".methodsNamespace", where, baseenv())
  error: cannot add binding of '.methodsNamespace' to the base environment
Fatal error: unable to initialize the JIT
```
