
# To build and push the docker container.

```
# from the root directory
docker build -f docker/g4v11.1.1-env -t dboogert/g4v11.1.1 .
docker push dboogert/g4v11.1.1:latest
```