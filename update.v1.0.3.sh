#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
echo $strPath

## update the biogen to VIP in app.py
sed -i "s|biogen|VIP|g" "$strPath/server/app/app.py"

## update the cellxgene title to cellxgene VIP
sed -i "s|cell&times;gene|cellxgene VIP|" "cellxgene/client/index_template.html"
sed -i "s|title=\"cellxgene|title=\"cellxgene VIP|" "cellxgene/client/src/components/app.js"
sed -i "s|  gene|  gene VIP|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"

sed -i "s|PLOTTING PANEL|Visualization in Plugin|i" "cellxgene/client/index_template.html"

cd cellxgene/client; make build; cp build/index.html $strPath/server/common/web/templates/

cd ../..

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strPath/server/common/web/static/.




