#!/usr/bin/env bash
## provide the location of python 3.7 lib install path as parameter
## e.g. config.sh
pythonV="$(python --version)"
if [[ $pythonV != *"Python 3.7"* && $pythonV != *"Python 3.8"* ]]; then
  echo "Only support Python 3.7 or 3.8"
  exit 0
fi

## obtain a clean version cellxgene a specific version by sha key
rm -rf cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene;git checkout 735eb11eb78b5e6c35ba84438970d0ce369604e1;cd ..

## update the source code for the VIP -----
echo -e "\nwindow.store = store;" >> cellxgene/client/src/reducers/index.js
read -d '' insertL << EOF
<link href='static/jspanel/dist/jspanel.css' rel='stylesheet'>
<!-- jsPanel JavaScript -->
<script src='static/jspanel/dist/jspanel.js'></script>
<!-- optional jsPanel extensions -->
<script src='static/jspanel/dist/extensions/modal/jspanel.modal.js'></script>
<script src='static/jspanel/dist/extensions/tooltip/jspanel.tooltip.js'></script>
<script src='static/jspanel/dist/extensions/hint/jspanel.hint.js'></script>
<script src='static/jspanel/dist/extensions/layout/jspanel.layout.js'></script>
<script src='static/jspanel/dist/extensions/contextmenu/jspanel.contextmenu.js'></script>
<script src='static/jspanel/dist/extensions/dock/jspanel.dock.js'></script>
<script>
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
        contentSize: {
            width: function() { return Math.min(730, window.innerWidth*0.9);},
            height: function() { return Math.min(400, window.innerHeight*0.5);}
        },
        position:    'center-top 0 180',
        animateIn:   'jsPanelFadeIn',
        headerControls:{
          close: 'remove',
          maximize: 'remove'
        },
        footerToolbar: '<span style="display:block; width:100%; height:4px; background-color:#F88519"></span>',
        headerTitle: 'Visualization in Plugin',
        contentAjax: {
            url: 'static/interface.html',
            done: function (panel) {
                   setInnerHTML(panel.content, this.responseText);
            }
        },
        onbeforeclose: function () {
            return confirm('Do you really want to close the panel?');
        }
    });
    plotPanel.minimize();
</script>
EOF
insertL=$(sed -e 's/[&\\/]/\\&/g; s/$/\\/' -e '$s/\\$//' <<<"$insertL")
sed -i "s|<div id=\"root\"></div>|$insertL\n&|g" "cellxgene/client/index_template.html" 

echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py
    
## update the cellxgene title to cellxgene VIP
sed -i "s|cell&times;gene|cellxgene VIP|" "cellxgene/client/index_template.html"
sed -i "s|title=\"cellxgene|title=\"cellxgene VIP|" "cellxgene/client/src/components/app.js"
sed -i "s|  gene|  gene VIP|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"


## buld the cellxgene and install -----------
conda remove PyYAML
conda install fsspec=0.6.3
cd cellxgene
make pydist
make install-dist
pip install 'scanpy==1.4.6'
cd ..

## finished setting up ------ 
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
strweb="${strPath}/server/common/web/static/."
echo $strweb
cp interface.html $strweb
cp jquery.min.js $strweb
cp -R DataTables $strweb
cp -R jspanel $strweb
cp cellxgene/server/test/decode_fbs.py $strPath/server/app/.
cp VIPInterface.py $strPath/server/app/.



