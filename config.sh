#!/usr/bin/env bash
## provide the location of python 3.7 lib install path as parameter
## e.g. config.sh
pythonV="$(python --version)"
if [[ $pythonV != *"Python 3.7"* && $pythonV != *"Python 3.8"* ]]; then
  echo "Only support Python 3.7 or 3.8"
  exit 0
fi

conda install -c conda-forge nodejs

## obtain a clean version cellxgene a specific version by sha key
rm -rf cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene;git checkout 735eb11eb78b5e6c35ba84438970d0ce369604e1;cd ..

## update the client-side source code of cellxgene for VIP
echo -e "\nwindow.store = store;" >> cellxgene/client/src/reducers/index.js
read -d '' insertL << EOF
<script src="static/jquery.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.6/d3.min.js" charset="UTF-8"></script>
<link href='static/jspanel/dist/jspanel.css' rel='stylesheet'>
<script src='static/jspanel/dist/jspanel.js'></script>
<script src='static/jspanel/dist/extensions/modal/jspanel.modal.js'></script>
<script src='static/jspanel/dist/extensions/tooltip/jspanel.tooltip.js'></script>
<script src='static/jspanel/dist/extensions/hint/jspanel.hint.js'></script>
<script src='static/jspanel/dist/extensions/layout/jspanel.layout.js'></script>
<script src='static/jspanel/dist/extensions/contextmenu/jspanel.contextmenu.js'></script>
<script src='static/jspanel/dist/extensions/dock/jspanel.dock.js'></script>
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
            url: 'static/interface.html',
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

sed -i "s|globals.datasetTitleMaxCharacterCount|50|; s|width: \"190px\"|width: \"300px\"|; s|{aboutURL ? <a href={aboutURL}|{myURL ? <a href={myURL}|; s|return|var myURL=displayTitle.split('_')[0].startsWith('GSE') \? 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='\+displayTitle.split('_')[0]:aboutURL;\n    \n    return|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"

sed -i "s|logoRelatedPadding = 50|logoRelatedPadding = 60|" "cellxgene/client/src/components/leftSidebar/index.js"

## update the cellxgene title to cellxgene VIP
sed -i "s|title=\"cellxgene\"|title=\"cellxgene VIP\"|" "cellxgene/client/src/components/app.js"

## update the server-side source code of cellxgene for VIP
echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py
    

## buld the cellxgene and install -----------
conda remove PyYAML
conda install fsspec=0.6.3
pip install tensorflow==2.2.0
pip install diffxpy==0.7.4

pip install plotly==4.8.1
pip install anndata==0.7.4
git clone https://github.com/theislab/scanpy.git
cd scanpy;git checkout 2ea9f836cec6e12a5cdd37bc4a229d4eadf59d37;cd ..
pip install scanpy/
pip install jupyter_client
pip install jupytext
pip install nbconvert
pip install pyarrow

# old versions
# pip install git+https://github.com/theislab/scanpy.git@split_show
# pip install 'scanpy==1.4.6'   # works for v1.4.6 too
# conda install -c plotly plotly-orca

cd cellxgene
make pydist
make install-dist
cd ..

## finished setting up ------ 
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
strweb="${strPath}/server/common/web/static/."
echo $strweb
cp interface.html $strweb
cp jquery.min.js $strweb
cp color_*.png $strweb

cp -R DataTables $strweb
cp -R jspanel $strweb

cp cellxgene/server/test/decode_fbs.py $strPath/server/app/.
cp VIPInterface.py $strPath/server/app/.

cp jquery-ui.min.js $strweb
cp color_*.png $strweb
cp -R ace $strweb
cp -R stackedbar $strweb
cp volcano.R $strPath/server/app/.
cp Density2D.R $strPath/server/app/.
