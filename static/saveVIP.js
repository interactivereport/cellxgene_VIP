/*
The following provide the functionality for saving and loading the session
*/

// functions for selected cells
function cellSave(){
  var v = {};
  for(const i of ['1','2']){
    var cN = window.store.getState().differential['celllist'+i];
    if(cN===null) cN = [];
    v['celllist'+i] = Array.from(cN);
  }
  return v;
}
function cellLoad(v){
  for(const i of [1,2]){
    window.store.dispatch({type: "store current cell selection as differential set "+i, data:Int32Array.from(v['celllist'+i])});
  }
}

// functions for genes and gene sets
function geneSave(){
  return {'genes':window.biogen,
          'geneSets':window.bioGeneGrp}
}
function geneLoad(v){
  window.biogen = v['genes'];
  window.bioGeneGrp = v['geneSets'];
}

// functions save all select/dropdown box
function selectSave(){
  var v = {};
  $('select').each(function(){
    if(!this.id.includes('geneset') && this.id.length>0){
      v[this.id] = $(this).val();
    }
  });
  return v;
}
function selectLoad(v){
  for(const one of Object.keys(v)){
    if(!!v[one]){
      $("#"+one).val(v[one]);
    }
  }
}

// function save all checkbox
function checkSave(){
  var cIDs = [];
  $('input:checkbox').each(function(){
    if(!!this.className && this.className.length>0){
      cIDs.push(this.className);
    }
  });
  var v={};
  for(const cID of [...new Set(cIDs)]){
    var notCK = [];
    for(const one of $("."+cID+":checkbox:not(:checked)")){
      notCK.push(one.value);
    }
    v[cID] = notCK;
  }
  return v;
}
function checkLoad(v){
  for(const cID of Object.keys(v)){
    $("."+cID).prop("checked",false);
    for(const one of $("."+cID)){
      if(!v[cID].includes(one.value)){
        if($(one).prop('disabled')){
          break;
        }else{
          $(one).trigger('click');
        }
      }
    }
  }
}

// function save all radio buttons
function radioSave(){
  var cNames = [];
  $('input:radio').each(function(){
    if(!!this.name && this.name.length>0){
      cNames.push(this.name);
    }
  });
  var v={};
  for(const cName of [...new Set(cNames)]){
    v[cName] = $("input[name='"+cName+"']:checked").val();
  }
  return v;
}
function radioLoad(v){
  for(const cName of Object.keys(v)){
    $("input[name='"+cName+"'][value='"+v[cName]+"']").prop('checked', true);
  }
}

// function save all input numbers
function numberSave(){
  var v={};
  $('input[type="number"').each(function(){
    v[this.id] = $(this).val();
  });
  return v;
}
function numberLoad(v){
  for(const cID of Object.keys(v)){
    $("#"+cID).val(v[cID]);
  }
}

// function save sortableUI
function sortableSave(){
  var v= {};
  $(".sortable-list").each(function(){
    var one = {}
    for(const eID of $(this).sortable('toArray')){
      one[eID]=$('[id="'+eID+'"]').find("label").text();
    }
    v[this.id] = one;
  });
  return v;
}
function sortableLoad(v){
  for(const cID in v){
    var htmlCheck = "";
    for(const one in v[cID]){
      htmlCheck += "<li id='"+one+"'><label>"+v[cID][one]+"</label></li>";
    }
    $("#"+cID).html(htmlCheck);
  }
}

// function save images
function imageCreate(){
  var v={};
  $('img').each(function(){
    if(this.id.length>0){
      v[this.id] = $(this).attr('src');
    }
  })
  return v;
}
function imageSave(){
  var v=imageCreate();

  var hiddenE = document.createElement('a');
  hiddenE.href="data:attachment/text,"+encodeURIComponent(JSON.stringify(v));
  hiddenE.target='_blank';
  hiddenE.download='cellxgene.'+window.store.getState().config.displayNames.dataset+'.img.txt';
  hiddenE.click();
}

// function save specific eID
function eIDsave(eIDs){
  var v={};
  for(const one of eIDs){
    switch(one){
      case 'GSEAenable':
      case "IMGcumuCK":
      case "GSEAcollapse":
        v[one]=$("#"+one).prop('checked');
        break;
    }
  }
  return v;
}
function eIDload(v){
  for(const one in v){
    switch(one){
      case 'GSEAenable':
      case "IMGcumuCK":
      case "GSEAcollapse":
        if(v[one]!=$("#"+one).prop('checked')){
          $("#"+one).trigger('click')
        }
        break;
    }
  }
}
