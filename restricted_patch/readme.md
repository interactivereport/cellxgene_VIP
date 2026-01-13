# Restricted patch

## Overview

This directory includes tools for installing the restricted patch.

This patch adds two new features to a **functional** VIP *conda* environment:
- allows disabling "Get Data" button and "Command Line Interface" tab
- allows adding identifier to timestamp files

Once installed, the patched *VIP* will have two additional flags available; however, it shall run just like before if these two flags are not specified for the instance.

## To install the patch

Execute:
```
$ bash install_VIPlight_restricted-patch.sh </path/to/VIP/environment>
```
for example:
```
$ bash install_VIPlight_restricted-patch.sh /anaconda/envs/VIP
```
You should see `Patching completed!` information once it's done.

## To remove the patch

Execute:
```
$ bash remove_restricted-patch.sh </path/to/VIP/environment>
```

You should see `Restricted Patch removed. Environment restored.` information once it's done.

If the environment was not patched previously, then nothing would be done.

## Usage

1. To run a VIP instance in restriction mode, simply add `--restricted` to the command. For example:
    ```
    $ /bin/cellxgene_VIP/VIPlight launch sample.h5ad --host 0.0.0.0 -p xxxx --max-category-items 500 --backed --disable-annotations --disable-gene-sets-save --restricted
    ```
    In the restriction mode, the "Get Data" button will always be disabled for Heatmap and Violin tabs, which prevents end user from downloading *csv* files. Additionally, the "Command Line Interface" will also be removed.

2. To run a VIP instance with identifier, simply provide a parameter `--identifier xyz` to the command. For example:
    ```
    $ /bin/cellxgene_VIP/VIPlight launch sample.h5ad --host 0.0.0.0 -p xxxx --max-category-items 500 --backed --disable-annotations --disable-gene-sets-save --identifier case1
    ```
    When running without identifier, there will be a `sample.timestamp` file created at the same directory of the `sample.h5ad` file for any user activities. In this example, the timestamp file will be `sample_case1.timestamp`. This would facilitate resource management, especially multi-instances on the same *h5ad*.

3. These two flags may be used together.
