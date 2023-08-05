# Building 3rd Party Dependencies

```
cd $G3FX_ROOT/3rdparty
mkdir .build
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$G3FX_ROOT/.install
make -j8
```