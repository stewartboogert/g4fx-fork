# Building
```
mkdir .build && cd .build 
cmake .. -DCMAKE_PREFIX_PATH=$PWD/../3rdparty/.install
```

# Testing todo

* output xunit xml
  * (catch2 output) https://github.com/catchorg/Catch2/blob/devel/docs/reporters.md
  * (circle read and output) https://circleci.com/docs/collect-test-data/
