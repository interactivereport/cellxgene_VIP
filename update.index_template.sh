#!/usr/bin/env bash
## obtain original index_template.html etc.
cd cellxgene
#git checkout bedbc87ed6178cd00a586feac3e99d4912d1c74e client/index_template.html
git checkout bdfd9fe0a5462a0c139675fe10356765d2bbd95b client/index_template.html
cd ..

sed -i "s|<div id=\"root\"></div>|$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' index_template.insert)\n&|" "cellxgene/client/index_template.html"
# The following line is for debug purpose
# sed -i "s|switch|console.log(action);\n  switch|" "cellxgene/client/src/reducers/graphSelection.js"

strPath=$(python -c "import server as _; print(_.__file__.replace('/server/__init__.py',''))")
cd cellxgene/client; make build
cp build/index.html $strPath/server/common/web/templates/
rm $strPath/server/common/web/static/main-*.*
rm $strPath/server/common/web/static/obsolete-*.*
cp build/static/*   $strPath/server/common/web/static/
cd ../..
