#!/usr/bin/env bash

# reading first position parameter as the VIP environment to remove patch
appEnvPath="${1:-~/.conda/envs/VIP_Du}"
appEnvPath=$(realpath ${appEnvPath/\~/$HOME})

echo -e "\nPatched location (first position parameter): $appEnvPath"
echo -e "Warning: changes will be made to the cellxgene of this environment!"
echo "Wait for 20 seconds, use Ctrl+C to terminate"
for i in {20..1}; do echo -ne "$i\033[0K\r"; sleep 1; done;

source $appEnvPath/etc/profile.d/conda.sh
conda activate
which python

strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"

cd $strPath/server

# list of files to restore
declare -A files=(
    ["default_config.py"]="default_config_prepatch.py"
    ["common/config/dataset_config.py"]="common/config/dataset_config_prepatch.py"
    ["common/config/server_config.py"]="common/config/server_config_prepatch.py"
    ["cli/launch.py"]="cli/launch_prepatch.py"
    ["common/web/templates/index.html"]="common/web/templates/index_prepatch.html"
)

# extra symlink
extra_symlink="common/web/static/interface_restricted.html"

# check if all prepatch files exist
missing=()
for orig in "${!files[@]}"; do
    pre="${files[$orig]}"
    if [ ! -e "$pre" ]; then
        missing+=("$pre")
    fi
done

if [ ${#missing[@]} -ne 0 ]; then
    echo "ERROR: cannot restore patch. Missing prepatch files:"
    for f in "${missing[@]}"; do
        echo "  - $f"
    done
    echo "Please verify your environment before proceeding."
    exit 1
fi

# restore each file
for orig in "${!files[@]}"; do
    pre="${files[$orig]}"
    echo "Restoring $orig from $pre ..."
    rm "$orig"
    mv "$pre" "$orig"
done

# remove the extra symlink
if [ -L "$extra_symlink" ]; then
    echo "Removing $extra_symlink..."
    rm "$extra_symlink"
fi

echo "Restricted Patch removed. Environment restored."
