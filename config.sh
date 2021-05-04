#!/usr/bin/env bash
## provide the location of python 3.7 lib install path as parameter
## e.g. config.sh


#envName="VIP"
#if [ $# -eq 0 ]
#then
#	# no user env name specified, use default
#    echo "No user env name provided, using default name VIP"
#else
#	envName=$1
#	echo "User provided env name $1"
#fi
#
## setup conda env based on VIP.yml
#
#conda env create -n $envName -f VIP.yml
#conda activate $envName
#echo "Done with conda env setup"
#

set -o errexit

pythonV="$(python --version)"
if [[ $pythonV != *"Python 3.7"* && $pythonV != *"Python 3.8"* ]]; then
  echo "Only support Python 3.7 or 3.8"
  exit 0
fi

## buld the cellxgene and install -----------

## obtain a clean version cellxgene a specific version by sha key
rm -rf cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene
git checkout bedbc87ed6178cd00a586feac3e99d4912d1c74e # v 0.16.7  # 735eb11eb78b5e6c35ba84438970d0ce369604e1 (v0.15.0)
sed -i 's|anndata>=0.7.0|anndata>=0.7.4|' 'server/requirements.txt'
sed -i 's|scanpy==1.4.6|scanpy==1.6.1|' 'server/requirements.txt'
cd ..

## update the client-side source code of cellxgene for VIP
echo -e "\nwindow.store = store;" >> cellxgene/client/src/reducers/index.js
read -d '' insertL << EOF
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jquery.min.js"></script>
<script src="https://d3js.org/d3.v4.min.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/stackedbar/d3.v3.min.js"></script>
<link href="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/jspanel.css" rel="stylesheet">
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/jspanel.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/modal/jspanel.modal.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/tooltip/jspanel.tooltip.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/hint/jspanel.hint.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/layout/jspanel.layout.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/contextmenu/jspanel.contextmenu.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/dock/jspanel.dock.js"></script>
<script>
    // execute JavaScript code in panel content
    var setInnerHTML = function(elm, html) {
        elm.innerHTML = html;
        Array.from(elm.querySelectorAll('script')).forEach( oldScript => {
            const newScript = document.createElement('script');
            Array.from(oldScript.attributes)
            .forEach( attr => newScript.setAttribute(attr.name, attr.value) );
            newScript.appendChild(document.createTextNode(oldScript.innerHTML));
            oldScript.parentNode.replaceChild(newScript, oldScript);
        });
    }
    var plotPanel = jsPanel.create({
        panelSize: '190 0',
        position: 'left-top 160 6',
        dragit: { containment: [-10, -2000, -4000, -2000] }, // set dragging range of VIP window
        boxShadow: 1,
        border: "solid #D4DBDE thin",
        contentOverflow: 'scroll scroll', // adding scrolling bars
        headerControls:{
          close: 'remove',
          minimize: 'remove',
          maximize: 'remove'
        },
        headerTitle: function () {return '<strong>Visualization in Plugin</strong>'},
        contentAjax: {
            url: window.location.href.replace(/\\\/+$/,'')+'/static/interface.html',
            done: function (panel) {
                   setInnerHTML(panel.content, this.responseText);
            }
        },
        onwindowresize: function(event, panel) {
            var jptop = parseInt(this.currentData.top);
            var jpleft = parseInt(this.currentData.left);
            if (jptop<-10 || window.innerHeight-jptop<10 || window.innerWidth-jpleft<10 || jpleft+parseInt(this.currentData.width)<10) {
                this.reposition("left-top 160 6");
            }
        },
        onunsmallified: function (panel, status) {
            this.reposition('center-top -370 180');
            this.resize({ width: 740, height: function() { return Math.min(480, window.innerHeight*0.6);} });
        },
        onsmallified: function (panel, status) {
            this.reposition('left-top 160 6');
            this.style.width = '190px';
        }
    }).smallify();
    plotPanel.headerbar.style.background = "#D4DBDE";
</script>
EOF
insertL=$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' <<<"$insertL")
sed -i "s|<div id=\"root\"></div>|$insertL\n&|" "cellxgene/client/index_template.html"

# sed -i "s|globals.datasetTitleMaxCharacterCount|50|; s|width: \"190px\"|width: \"300px\"|; s|{aboutURL ? <a href={aboutURL}|{myURL ? <a href={myURL}|; s|return|var myURL=displayTitle.split('_')[0].startsWith('GSE') \? 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='\+displayTitle.split('_')[0]:aboutURL;\n    \n    return|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"

sed -i "s|logoRelatedPadding = 50|logoRelatedPadding = 60|" "cellxgene/client/src/components/leftSidebar/index.js"

## update the cellxgene title to cellxgene VIP
sed -i "s|title=\"cellxgene\"|title=\"cellxgene VIP\"|" "cellxgene/client/src/components/app.js"

## modify zoom/pan default
sed -i "s|const *scaleMax *= *[0-9\.]\+|const scaleMax = 50000|; s|const *scaleMin *= *[0-9\.]\+|const scaleMin = 0.1|; s|const *panBound *= *[0-9\.]\+|const panBound = 80|" "cellxgene/client/src/util/camera.js"

## update the server-side source code of cellxgene for VIP
## Please change /tmp if different temporary directory is used in your local environment
echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py


# Old branch for nicer plots that are incorporated into ver 1.6.1 now
#git clone https://github.com/theislab/scanpy.git
#cd scanpy;git checkout 2ea9f836cec6e12a5cdd37bc4a229d4eadf59d37;cd ..
#pip install scanpy/

cd cellxgene
make pydist
make install-dist
cd ..

## finished setting up ------
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb
cp volcano.R $strPath/server/app/.
cp violin.R $strPath/server/app/.
cp Density2D.R $strPath/server/app/.
cp bubbleMap.R $strPath/server/app/.
cp fgsea.R $strPath/server/app/.
mkdir -p $strPath/server/app/gsea
cp gsea/*gmt $strPath/server/app/gsea
cp vip.env $strPath/server/app/. 2>/dev/null | true
sed -i "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;

echo -e "\nls -l $strweb\n"
ls -l $strweb

export LIBARROW_MINIMAL=false
if [ $(python -c 'import nbconvert; print(nbconvert.__version__)') != "5.6.1" ]; then pip install nbconvert==5.6.1; fi
