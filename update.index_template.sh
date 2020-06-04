#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
echo $strPath

## obtain original index_template.html

cd cellxgene; git checkout 735eb11eb78b5e6c35ba84438970d0ce369604e1 client/index_template.html client/src/components/leftSidebar/topLeftLogoAndTitle.js client/src/components/leftSidebar/index.js; cd ..

## update the source code for the VIP -----
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
        panelSize: '190 0',
        position:    'left-top 160 6',
        animateIn:   'jsPanelFadeIn',
        boxShadow:    1,
        border: "solid #AFBEC4 thin",
        contentOverflow: 'scroll scroll',
        headerControls:{
          close: 'remove',
          minimize: 'remove',
          maximize: 'remove'
        },
        footerToolbar: '<span style="display:block; width:100%; height:4px; background-color:#AFBEC4"></span>',
        headerTitle: 'Visualization in Plugin',
        contentAjax: {
            url: 'static/interface.html',
            done: function (panel) {
                   setInnerHTML(panel.content, this.responseText);
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
</script>
EOF
insertL=$(sed -e 's/[&\\/]/\\&/g; s/$/\\/' -e '$s/\\$//' <<<"$insertL")
sed -i "s|<div id=\"root\"></div>|$insertL\n&|;s|cell&times;gene|cellxgene VIP|" "cellxgene/client/index_template.html"

sed -i "s|globals.datasetTitleMaxCharacterCount|50|; s|width: \"190px\"|width: \"300px\"|; s|{aboutURL ? <a href={aboutURL}|{myURL ? <a href={myURL}|; s|return|var myURL=displayTitle.split('_')[0].startsWith('GSE') \? 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='\+displayTitle.split('_')[0]:aboutURL;\n    \n    return|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"

sed -i "s|logoRelatedPadding = 50|logoRelatedPadding = 60|" "cellxgene/client/src/components/leftSidebar/index.js"

cd cellxgene/client; make build; cp build/index.html $strPath/server/common/web/templates/; cd ../..
