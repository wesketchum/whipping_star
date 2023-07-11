podman machine stop podman-machine-default
podman machine rm podman-machine-default
podman machine init -v $PWD:$PWD
podman machine set --memory 8192
podman machine start
