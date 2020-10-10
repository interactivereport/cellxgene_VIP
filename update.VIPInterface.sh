#!/usr/bin/env bash
if [ -n "$1" ]; then
echo "usually update once"
fi

## finished setting up ------
strPath=$(python -c "import server as _; print(_.__file__.replace('/server/__init__.py',''))")
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb

if [ -n "$1" ]; then
cp -R stackedbar $strweb
cp -R d3plot $strweb
find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;
fi

echo -e "\nls -l $strweb\n"
ls -l $strweb
