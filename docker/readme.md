# To build and push the docker container.


Assumes $G4FX_ROOT is the clone location f
```
cd $G3FX_ROOT
docker build -f docker/g4v11.1.1-env -t dboogert/g4v11.1.1 .
docker push dboogert/g4v11.1.1:latest
```