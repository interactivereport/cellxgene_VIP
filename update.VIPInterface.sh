#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
strweb="${strPath}/server/common/web/static/."
echo $strPath

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strPath/server/common/web/static/.
cp -R stackedbar $strPath/server/common/web/static/.

if [ -n "$1" ]; then
	echo "usually update once"
	cp jquery-ui.min.js $strweb
	cp color_*.png $strweb
	cp -R ace $strweb
	cp -R stackedbar $strweb
	cp volcano.R $strPath/server/app/.
	cp Density2D.R $strPath/server/app/. 
	pip install plotly==4.8.1
	pip install jupytext
	pip install nbconvert
fi
