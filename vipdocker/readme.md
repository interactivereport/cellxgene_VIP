## Build Image

ensure you have docker access
```
$ command -v "docker"
        /usr/bin/docker
$ docker ps
        CONTAINER ID   IMAGE                           COMMAND                  CREATED        STATUS        PORTS
```

navigate to repo main page (/path/to/BxGenomics/scRNAview or /path/to/cellxgene_VIP)
```
$ cd </path/to/>cellxgene_VIP
```
execute building
```
$ nohup docker buildx build -t bxgenomics_vip . > vipdocker/build/build.log 2>&1 &
```
check the content of `vipdocker/build/build.log` to monitor the building progress

This will create a docker image named "bxgenomics_vip"


## Prepare to run

ensure that the image is available locally:
```
$ docker images | grep "bxgenomics_vip"
```

if unavailable locally, build it or pull it from docker hub:
```
$ docker pull dujiang1031/bxgenomics_vip
```
add the wrapper to PATH:
```
$ export PATH=</path/to/cellxgene_VIP>/vipdocker/run:$PATH
```

## Run

navigate to your working directory, first execute
```
$ VIPdocker </path/to/working/dir>
```
this will generate a **run.yml** file. Open and edit it:
```
vip_docker_name: bxgenomics_vip
docker_host_mount:
  - /share
  - /aws_s3

app_uid: "1014:1999" # userid:groupid; make sure that this user has write permission on the folder of sample.h5ad
tmp_location: /path/to/working/dir # location to create a /tmp folder for cellxgene

h5ad: /absolute/path/to/sample.h5ad
port: "8212"

# default args already included in the docker and not needed here:
#   --host=0.0.0.0
#   --max-category-items=500
#   --backed
#   --disable-annotations
#   --disable-gene-sets-save
# additional args can be provided here. Note that "=" is required for parameters
launch_args:
  - --restricted
  - --identifier=xyz

```
after modifying and saving the **run.yml** file, run it by:
```
$ VIPdocker run.yml
```
or
```
$ nohup VIPdocker run.yml > run.log 2>&1 &
```
then navigate to localhost:<\port> for visualization & analysis
