The following LEfSe error is caused by newer versions of `rpy2`.

```
AttributeError: 'NoneType' object has no attribute 'rownames'
```

`rpy2` needs to be downgraded to `3.5.10`.

```bash
pip install rpy2==3.5.10
```
