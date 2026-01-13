#!/usr/bin/env bash

# reading first position parameter as the VIP environment to be patched
appEnvPath="${1:-~/.conda/envs/VIP_Du}"
appEnvPath=$(realpath ${appEnvPath/\~/$HOME})

echo -e "\nPatching location (first position parameter): $appEnvPath"
echo -e "Warning: changes will be made to the cellxgene of this environment!"
echo "Wait for 20 seconds, use Ctrl+C to terminate"
for i in {20..1}; do echo -ne "$i\033[0K\r"; sleep 1; done;

patchPath=$(readlink -e $(dirname $0))
echo $patchPath

source $appEnvPath/etc/profile.d/conda.sh
conda activate
which python

strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"

cd $strPath/server
ln -s $patchPath/interface_restricted.html common/web/static/interface_restricted.html
mv default_config.py default_config_prepatch.py
ln -s $patchPath/default_config.py default_config.py
mv common/config/dataset_config.py common/config/dataset_config_prepatch.py
ln -s $patchPath/dataset_config.py common/config/dataset_config.py
mv common/config/server_config.py common/config/server_config_prepatch.py
ln -s $patchPath/server_config.py common/config/server_config.py
mv cli/launch.py cli/launch_prepatch.py
ln -s $patchPath/launch.py cli/launch.py
mv app/app.py app/app_prepatch.py
ln -s $patchPath/app.py app/app.py
mv common/web/templates/index.html common/web/templates/index_prepatch.html
ln -s $patchPath/index.html common/web/templates/index.html

echo "Patching completed!"
