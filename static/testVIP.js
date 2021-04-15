/*
The following provide the functionality for testing the VIP after updating using MS data
*/
var testVIPimg, testVIPlist;

function testVIP(eID){
  // if the test has started already
  $(".tVIPbt").prop('disabled',true);
  $("#"+eID).html("");
  testVIPhtmladd(eID,"<h3>TESTING VIP</h3>This might take some time, please be patient.<br>");

  testVIPgetData('static/testVIP/'+window.store.getState().config.displayNames.dataset+'.info.txt',testVIPload,eID);
  testVIPgetData('static/testVIP/'+window.store.getState().config.displayNames.dataset+'.img.txt',testVIPimgGet,eID);
  //All plots testing sequentially in testVIPimgGet function, since the async data retriving
}

function testVIPhtmladd(eID,message){
  $("#"+eID).html($("#"+eID).html()+message);
}
function testVIPgetData(url,callback){
  var args = Array.prototype.slice.call(arguments, 2);
  var xhr = new XMLHttpRequest();
  xhr.onload = function(){
    if(xhr.readyState === 4 && xhr.status === 200){
      callback.apply(xhr,args);
    }
  }
  xhr.open('GET', url, true);
  xhr.send(null);
}
function testVIPload(eID){
  var content = JSON.parse(decodeURIComponent(this.responseText));
  loadexe(content);
  testVIPhtmladd(eID,"Loading setting complete<br>");
}
function testVIPimgGet(eID){
  testVIPimg = JSON.parse(decodeURIComponent(this.responseText));
  testVIPhtmladd(eID,"Loading information complete<br>");
  testVIPlist = Object.keys(testVIPimg);
  testVIPnext(eID);
}
// The following testing all of them at once cause server 'Segmentation fault'
// Segmentation fault      (core dumped) cellxgene launch ...
function testVIPall(eID){
  for(const one in testVIPimg){
    var fn = window[one.replace('img','plot')];
    if(one.includes("DEG")){
      fn = window["DEGfind"];
    }
    if(typeof fn === "function") fn();
    testVIPhtmladd(eID,"Testing"+one.replace('img','')+" ...<br>");
    testVIPone(one,eID,0);
  }
}
function testVIPnext(eID){
  if(testVIPlist.length>0){
    var one = testVIPlist.pop();
    var fn = window[one.replace('img','plot')];
    if(one=="DEGimg"){
      fn = window["DEGfind"];
    }else if(one=="preDEGvolcanoimg"){
      fn = window["preDEGvolcanofind"];
    }
    if(typeof fn === "function"){
      $("#"+one).remove();
      testVIPhtmladd(eID,"Testing "+one.replace('img','')+" ...<br>");
      if(one=="preDEGVolcanoimg"){
        fn('preDEGvolcano');
      }else{
        fn();
      }
      testVIPone(one,eID,0);
    }else{
      testVIPnext(eID);
    }
  }else{
    testVIPhtmladd(eID,"=========== Testing complete! ===========");
    $(".tVIPbt").prop("disabled",false);
  }
}
function testVIPone(imgID,eID,nTimes){
  if(nTimes>60){//1 minutes
    testVIPhtmladd(eID,"<p style='color:red;'>Testing on "+imgID+" timeout</p>");
    testVIPnext(eID);
    return;
  }
  var info = $('#'+imgID).attr('src');
  if(info===undefined){
    setTimeout(testVIPone,1000,imgID,eID,nTimes+1);
  }else{
    if(info===testVIPimg[imgID]){
      testVIPhtmladd(eID,"Testing on "+imgID+" complete successfully<br>");
    }else if(info.length>100){
      testVIPhtmladd(eID,"<p style='color:salmon;'>Testing on "+imgID+" image changed</p>");
    }
    else{
      testVIPhtmladd(eID,"<p style='color:darkred;'>Testing on "+imgID+" failed</p>");
    }
    testVIPnext(eID);
  }
}

// the following are creating a test case
function createTest(grpName='',eID){
  if(grpName.length==0 || !Object.keys(window.store.getState().categoricalSelection).includes(grpName)){
    annoNames = Object.keys(window.store.getState().categoricalSelection);
    grpName = randomSel(annoNames,1)[0];
  }
  testVIPhtmladd(eID,"Creating test case on "+grpName+", please be patient, this might take a while. A successful message will be shown at the end.<br>");
  randomSelCell(grpName,eID);
  randomInitDEG(grpName,eID);
}

function randomSel(x,num){
  if(x.length<num){
    return x;
  }
  var ix = new Array(x.length).fill().map((a, i) => a = i).sort(() => Math.random() - 0.5).splice(0,num);
  var sel=ix.map(i=>x[i]);
  return sel;
}
function randomSelCell(grpName,eID){
  testVIPhtmladd(eID,"Randomly selecting cells ...<br>");
  var gUnique = window.store.getState().annoMatrix.schema.annotations.obsByName[grpName].categories;
  var gNames = randomSel(gUnique,~~(gUnique.length/2));
  gNames = [gNames,gUnique.filter(i=>!gNames.includes(i))];

  var gList = window.store.getState().annoMatrix._cache.obs.__columns[window.store.getState().annoMatrix._cache.obs.colIndex.getOffset(grpName)];
  var setN=1;
  for(const one of gNames){
    var ix=gList.map((e,i) => one.includes(e)?i:'').filter(String);
    window.store.dispatch({type: "store current cell selection as differential set "+setN, data:Int32Array.from(randomSel(ix,Math.min(500,~~(ix.length/10))))});
    setN++;
  }
}
function randomInitDEG(grpName,eID){
  testVIPhtmladd(eID,"Finding DEGs on selected cells ...<br>");
  DEGfind();
  setTimeout(randomSelGene,1000,grpName,eID);
}
function randomSelGene(grpName,eID){
  testVIPhtmladd(eID,". ");
  if(window.DEGraw === undefined){
    setTimeout(randomSelGene,1000,grpName,eID);
  }else{
    var tableD = d3.csvParse(window.DEGraw,function(data){
      return [].concat([data.gene,
                (+data.log2fc).toFixed(2),
                (+data.pval).toExponential(2),
                (+data.qval).toExponential(2)],
                Object.values(data).slice(5));
      });
    var ipos=0;
    window.biogen=[];
    for(var i=0;i<5;i++){
      var gName = [];
      for(var j=0;j<~~(Math.random()*10)+1;j++){
        gName.push(tableD[ipos+j][0]);
        window.biogen.push(tableD[ipos+j][0]);
      }
      ipos += gName.length;
      window.bioGeneGrp['set'+(1+i)]=gName;
    }
    sync();
    randomSelOpt(grpName,eID)
  }
}
function randomSelOpt(grpName,eID){
  randomSelBox(grpName,eID);
  randomSelDropdown(grpName,eID);
  setTimeout(randomPlot,500,eID);
}
function randomSelBox(grpName,eID){
  var cIDs = [];
  $('input:checkbox').each(function(){
    if(!!this.className && this.className.length>0){
      cIDs.push(this.className);
    }
  });
  testVIPhtmladd(eID,"<br>Randomize checkbox group: ");
  for(const cID of [...new Set(cIDs)]){
    if(cID.includes('ABBR') || cID.includes('CLI') || cID.includes('cell')) continue;
    testVIPhtmladd(eID,cID+"; ");
    $("."+cID).each(function(){
      if($(this).prop('checked')){
        $(this).trigger('click');
      }
    });
    for(const one of randomSel($("."+cID),Math.max(2,~~(Math.random()*$("."+cID).length/2)))){
      if($(one).prop('disabled')){
        break;
      }else{
        $(one).trigger('click');
      }
    }
  }
}
function randomSelDropdown(grpName,eID){
  testVIPhtmladd(eID,"<br>Randomize dropdown: ");
  $('select').each(function(){
    if(!this.id.includes('geneset') && this.id.length>0){
      testVIPhtmladd(eID,this.id+"; ");
      var opts=[],preSel=grpName;
      $("#"+this.id+" option").each(function() {
        opts.push($(this).val());
        if(/^Welch/.test($(this).val()) && /cellxgene$/.test($(this).val())){
          preSel = $(this).val();
        }
      });
      if(this.id.includes('split') && opts.includes('None')){
        preSel = 'None';
      }
      if(opts.includes(preSel)){
        $(this).val(preSel);
      }else{
        $(this).val(randomSel(opts,1));
      }
    }
  });
  sync();
}

function randomPlotOne(imgID,strBT,nTimes,eID,bSaved){
  if(nTimes>=120){
    testVIPhtmladd(eID,"<p style='color:darkred;'>Error: Timeout!</p>");
    $(".tVIPbtA").prop("disabled",false);
    enableVIPtest();
    return;
  }
  if(imgID.length>0 && $('#'+imgID).html().length<100){
    setTimeout(randomPlotOne,500,imgID,strBT,nTimes+1,eID,bSaved);
  }else{
    if(strBT.length==0){
      testVIPhtmladd(eID,"<br>Plotting is completed<br>");
      if(bSaved){
        setTimeout(randomSave,500,eID);
      }
    }else{
      var one = strBT.pop();
      var fn = window[one], imgID="";
      if(!one.includes('preDEG')){
        $("."+one).each(function(){
          if($(this).text().includes("Plot")){
            imgID=one.replace('bt','resize');
            var pID = $("#"+imgID).find('p').prop('id');
            if(pID===undefined){
              $("#"+imgID).html('');
            }else{
              $("#"+imgID).html('<p id="'+pID+'"></p>');
            }
            testVIPhtmladd(eID,one+"; ");
            $(this).trigger("click");
            return;
          }
        });
      }
      setTimeout(randomPlotOne,500,imgID,strBT,0,eID,bSaved);
    }
  }
}
function randomPlot(eID,bSaved=true){
  testVIPhtmladd(eID,"<br>Starting plotting<br>Executing: ");
  var strBT=[];
  $("button[class$='bt']").each(function(){
    strBT.push(this.className);
  });
  strBT = [...new Set(strBT)];
  setTimeout(randomPlotOne,500,"",strBT,0,eID,bSaved);
}
function randomSave(eID){
  testVIPhtmladd(eID,"<br>Saving ...<br>");
  var D={'method':'saveTest',
         'dataset': window.store.getState().config.displayNames.dataset+'.h5ad',
         'info':encodeURIComponent(JSON.stringify(saveContent())),
         'img':encodeURIComponent(JSON.stringify(imageCreate()))};
  $.ajax({
    type:"POST",
    url: VIPurl,
    data:JSON.stringify(D),
    contentType: 'application/json;charset=UTF-8',//
    success: function(res){
      testVIPhtmladd(eID,"Test information was saved successfully!");
      $(".tVIPbtA").prop("disabled",false);
      enableVIPtest();
    },
    error: function(request,status,error){
      testVIPhtmladd(eID,request.status+request.responseText);
      $(".tVIPbtA").prop("disabled",false);
      enableVIPtest();
    }
  });
}


function randomViolin(grpName){
  if($('#'+imgID).attr('src')===undefined){
    if($("#SGVspin").css('visibility')==="hidden"){
      $("#SGVgene").val(randomSel(window.biogen,1)[0]);
      $("#SGVgrp").val(grpName);
      SGVplot('plot');
    }
    setTimeout(randomViolin,500,grpName);
  }else{

  }
}
