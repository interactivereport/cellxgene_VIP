/*
The following provide the functionality for testing the VIP after updating using MS data
*/
var testVIPimg, testVIPlist;

function testVIP(eID){
  // if the test has started already
  if($("#"+eID).html().includes("TESTING VIP")){
    return;
  }
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
    $("#"+one).remove();
    testVIPhtmladd(eID,"Testing "+one.replace('img','')+" ...<br>");
    var fn = window[one.replace('img','plot')];
    if(one.includes("DEG")){
      fn = window["DEGfind"];
    }
    if(typeof fn === "function") fn();
    testVIPone(one,eID,0);
  }else{
    testVIPhtmladd(eID,"=========== Testing complete! ===========<br>");
    testVIPhtmladd(eID,"<button onclick='$(\"#"+eID+"\").html(\"\");testVIP(\""+eID+"\")'>Retest</button>");
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

function createTest(grpName=''){
  if(grpName.length==0 || Object.keys(window.store.getState().categoricalSelection).includes(grpName)){
    grpName = randomSel(Object.keys(window.store.getState().categoricalSelection),1)[0];
  }
  console.log("Creating test case, please be patient, this might take a while. A successful message will be shown at the end.");
  randomSelCell(grpName);
  DEGfind();
  setTimeout(randomSelGene,1000,grpName);
}

function randomSel(x,num){
  if(x.length<num){
    return x;
  }
  var ix = new Array(x.length).fill().map((a, i) => a = i).sort(() => Math.random() - 0.5).splice(0,num);
  var sel=ix.map(i=>x[i]);
  return sel;
}
function randomSelCell(grpName){
  console.log("Randomly selecting cells ...");
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
function randomSelGene(grpName){
  console.log("Selecting DEGs ...");
  if(window.DEGraw === undefined){
    setTimeout(randomSelGene,1000,grpName);
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
    randomSelOpt(grpName)
  }
}
function randomSelOpt(grpName){
  randomSelBox(grpName);
  randomSelDropdown(grpName);
  setTimeout(randomPlot(),500);
}
function randomSelBox(grpName){
  var cIDs = [];
  $('input:checkbox').each(function(){
    if(!!this.className && this.className.length>0){
      cIDs.push(this.className);
    }
  });
  for(const cID of [...new Set(cIDs)]){
    if(cID.includes('ABBR') || cID.includes('CLI') || cID.includes('cell')) continue;
    console.log("Randomize checkbox group "+cID);
    $("."+cID).each(function(){
      $(this).prop('checked', false);
    });
    for(const one of randomSel($("."+cID),Math.max(2,~~(Math.random()*$(".DENS2Dgenes").length/2)))){
      if($(one).prop('disabled')){
        break;
      }else{
        $(one).trigger('click');
      }
    }
  }
}
function randomSelDropdown(grpName){
  $('select').each(function(){
    if(!this.id.includes('geneset') && this.id.length>0){
      console.log("Randomize dropdown "+this.id);
      var opts=[],strDEGmethod="Welch's t-test (cellxgene)";
      $("#"+this.id+" option").each(function() { opts.push($(this).val()); });
      if(opts.includes(grpName)){
        $(this).val(grpName);
      }else if(opts.includes(strDEGmethod)){
        $(this).val(strDEGmethod);
      }else{
        $(this).val(randomSel(opts,1));
      }
    }
  });
  sync();
}

function randomPlotOne(imgID,strBT,nTimes){
  if(nTimes>=120){
    console.log("Error: Timeout!");
    return;
  }
  if(imgID.length>0 && $('#'+imgID).html().length<100){
    setTimeout(randomPlotOne,500,imgID,strBT,nTimes+1);
  }else{
    if(strBT.length==0){
      setTimeout(randomSave,500);
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
            console.log("Executing "+one);
            $(this).trigger("click");
            return;
          }
        });
      }
      setTimeout(randomPlotOne,500,imgID,strBT,0);
    }
  }
}
function randomPlot(){
  console.log("Starting plotting ...");
  var strBT=[];
  $("button[class$='bt']").each(function(){
    strBT.push(this.className);
  });
  strBT = [...new Set(strBT)];
  setTimeout(randomPlotOne,500,"",strBT,0);
}
function randomSave(){
  console.log("Saving ...");
  var D={'method':'saveTest',
         'info':encodeURIComponent(JSON.stringify(saveContent())),
         'img':encodeURIComponent(JSON.stringify(imageCreate()))};
  $.ajax({
    type:"POST",
    url: VIPurl,
    data:JSON.stringify(D),
    contentType: 'application/json;charset=UTF-8',//
    success: function(res){
      console.log("Test information was saved successfully!");
    },
    error: function(request,status,error){
      console.log(request.status+request.responseText);
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
