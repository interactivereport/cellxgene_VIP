#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
echo $strPath
cp biogenInterface.py $strPath/server/app/.
cp interface.html $strPath/server/common/web/static/.

read -d '' strScanpy << EOF
    _utils.savefig_or_show('dotplot', show=show, save=save)

EOF
read -d '' newScanpy << EOF
    _utils.savefig_or_show('dotplot', show=show, save=save)
    return fig
EOF
strScanpy=$(sed -e 's/[&\\/]/\\&/g; s/$/\\/' -e '$s/\\$//' <<<"$strScanpy")
newScanpy=$(sed -e 's/[&\\/]/\\&/g; s/$/\\/' -e '$s/\\$//' <<<"$newScanpy")

sed -i "s|$strScanpy|$newScanpy|g" "${strPath}/scanpy/plotting/_anndata.py"
