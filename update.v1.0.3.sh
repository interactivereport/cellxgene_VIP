#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
echo $strPath

## update the biogen to VIP in app.py
sed -i "s|biogen|VIP|g" "$strPath/server/app/app.py"

## update the cellxgene title to cellxgene VIP
sed -i "s|gene|geneVIP|g" "cellxgene/client/index_template.html"
sed -i "s|gene|geneVIP|g" "cellxgene/client/index.html"
sed -i "s|PLOTTING PANEL|Visualization in Plugin|g" "cellxgene/client/index_template.html"
cd cellxgene/client; make build; cd ..;make install-dist;cd ..;

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strPath/server/common/web/static/.




