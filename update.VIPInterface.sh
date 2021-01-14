#!/usr/bin/env bash
if [ -n "$1" ]; then
echo "usually update once"
fi

## finished setting up ------
strPath=$(python -c "import server as _; print(_.__file__.replace('/server/__init__.py',''))")
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb
cp vip.env $strPath/server/app/. 2>/dev/null | true

if [ -n "$1" ]; then
cp saveVIP.js $strweb
cp -R stackedbar $strweb
cp -R d3plot $strweb
cp bubbleMap.R $strPath/server/app/.
cp violin.R $strPath/server/app/.
sed -i "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;
fi

echo -e "\nls -l $strweb\n"
ls -l $strweb
