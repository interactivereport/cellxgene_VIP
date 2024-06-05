#!/usr/bin/env python

import sys, io, random, time, getopt, math
import anndata as ad
import pandas as pd
import numpy as np
#from base64 import b64encode, b64decode
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline

def findGenes(genes,gNames):
    genes = genes.lower().split(",")
    allG = {v.lower():v for v in gNames}
    return [allG[x] for x in sorted(set(allG.keys()) & set(genes),key=genes.index)]
def randomCell(selN,cNames):
    selC = [False]*len(cNames)
    for i in random.sample(list(range(len(cNames))),selN):
        selC[i] = True
    return selC
def checkGene(selG,genes):
    if len(selG)==0:
        print("ERROR: none of provided genes ("+",".join(genes if isinstance(genes,list) else [genes])+") listed in the h5ad file!")
        exit()
def checkGrp(grp,allGrp):
    if not grp in allGrp:
        print("ERROR: annotation ("+grp+") is unknown in the h5ad file!")
        exit()
def filterCell(df,gCut,logFlag=False):
    if gCut is not None:
        df = df[df[df>float(gCut)].count(axis=1)>0]
    if logFlag:
        df=np.log2(1+df)
    return df
def obtainData(strH5ad,genes,grp,gN=None,cN=None,gCut=None,logMax=float("inf"),vname=None):
    sT = time.time()
    D=ad.read_h5ad(strH5ad,backed='r',as_sparse=["X"],chunk_size=60000)
    if vname is not None and vname in D.var.columns:
      D.var_names = list(D.var[vname])
    checkGrp(grp,D.obs.columns)
    logFlag = False
    if logMax!=float("inf") and np.max(D.X.value)>logMax:
        logFlag = True
    selG = findGenes(genes,D.var_names)
    checkGene(selG,genes)
    if gN is not None:
        selG=selG[0:gN]
    if cN is None:
        X = pd.concat([filterCell(D[:,selG].to_df(),gCut,logFlag),D.obs[grp]],axis=1,join='inner')
    elif int(cN)<D.shape[0] and int(cN)>10:
        selC=randomCell(int(cN),D.obs_names)
        subD = filterCell(D[:,selG].to_df().loc[selC,:],gCut,logFlag)
        X = pd.concat([subD,D.obs.loc[selC,grp]],axis=1,join='inner')
    else:
        raise Exception("ERROR: input cell number ("+cN+") is invalid!")
    return X,selG

def violin(strH5ad,genes,grp,cN=None,gCut=None,vname=None):
    cG=1
    X,selG = obtainData(strH5ad,genes,grp,cG,cN=cN,gCut=gCut,vname=vname)
    selG=selG[0]
    fig = go.Figure() #layout=go.Layout(paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
    for one in X[grp].cat.categories:
        fig.add_trace(go.Violin(x=X[grp][X[grp]==one],
                                y=X[selG][X[grp]==one],
                                name=one,
                                box={"visible":True,"width":0.1,"fillcolor":'#ffffff',"line":{"color":'#000000'}},
                                meanline_visible=True,points=False))
    annoMaxL=max([len(str(i)) for i in X[grp].cat.categories])
    xlabH=annoMaxL*8
    xlabW=annoMaxL*7
    h=50+30*cG+xlabH
    w=160+30*len(X[grp].cat.categories)+xlabW
    fig.update_layout(plot_bgcolor='#fff',title={"text":grp},title_x=0.5,showlegend=False,
                        colorway=px.colors.qualitative.Set1,
                        height=h,width=w,
                        margin={"l":0,"r":0,"t":30,"b":0})
    fig.update_yaxes(range=[0,X[selG].max()],gridcolor="#eee",linecolor="#000",title={"text":selG,"font":{"color":"#000"}})
    fig.update_xaxes(linecolor="#000",tickfont={"size":15})
    #fig.write_html("test1.html")
    #fig = px.violin(X,y=selG[0],color=grp,box={"visible":True,"fillcolor":"ffffff","line":{"color":"000000"}},x=grp,range_y=[0,X[selG[0]].max()])
    buffer = io.StringIO()
    fig.write_html(buffer,include_plotlyjs=False,full_html=False)
    #html=buffer.getvalue()
    #divIndex=html.find('<div id')
    #if html.find('<div id',divIndex+1)!=-1:
    #    print("Error: plotly structure change! Contact oyoung@bioinforx.com")
    #    exit()
    #scriptIndex=html.find('</script>',divIndex)
    #imgD = buffer.getvalue()
    #imgD = b64encode(imgD.encode()).decode()
    #with open("test1.txt","w") as f:
    #    f.write(imgD)
    return buffer.getvalue() #html[divIndex:(scriptIndex+len("</script>"))]

def dot(strH5ad,genes,grp,cN=None,gCut=None,expRange=[None,None],permax=None,logmax=float("inf"),legendDirection="v",vname=None):
    X,selG = obtainData(strH5ad,genes,grp,cN=cN,gCut=gCut,logMax=logmax,vname=vname)
    df = list()
    for anno in X[grp].cat.categories:
        for gene in selG:
            tmp = X.loc[X[grp]==anno,gene]
            df += [[anno,gene,
                    round(tmp.median(),2),
                    round(tmp.mean(),2),
                    round(tmp[tmp>0].count()/tmp.shape[0]*100,2),
                    tmp.shape[0]]]
    DOT=pd.DataFrame(df,columns=['grp','gene','median','mean','percent','cellN'])
    DOT=DOT.fillna(0)
    percentSizeLegendW=80
    percentSizeLegendH=80
    annoMaxL=max([len(str(i)) for i in X[grp].cat.categories])
    xlabH=35+annoMaxL*5
    xlabW=annoMaxL*7
    h=50+30*len(selG)+xlabH
    w=160+30*len(X[grp].cat.categories)+xlabW
    if legendDirection=="v":
        w += percentSizeLegendW
        subF = {"dot":[1,1],"legend":[1,2]}
        wR = [1-percentSizeLegendW/w, percentSizeLegendW/w]
        hR = None
        legendXY={"x":[1,1,1],"y":[1,2,3]}
    else:
        h += percentSizeLegendH
        subF = {"dot":[2,1],"legend":[1,1]}
        wR = None
        hR = [percentSizeLegendH/h, 1-percentSizeLegendH/h]
        legendXY={"y":[1,1,1],"x":[1,2,3]}
    if permax is None:
        permax = (int)(DOT["percent"].max())
    figDot = make_subplots(rows=max(subF['dot'][0],subF['legend'][0]),
                           cols=max(subF['dot'][1],subF['legend'][1]),
                           column_widths=wR,row_heights=hR,
                           horizontal_spacing=0,vertical_spacing=0,
                           subplot_titles=("","Percentage"))
    figDot.update_annotations(font_size=11)
    figDot.add_trace(go.Scatter(x=DOT["grp"],y=DOT["gene"],mode='markers',
            text=['%s'%i for i in DOT['cellN']],
            name="",
            marker=dict(size=DOT['percent'],color=DOT["mean"],
                cmin=expRange[0],cmax=expRange[1],
                colorscale=['rgb(0,0,255)','rgb(150,0,90)','rgb(255,0,0)'],
                sizemin=0,
                sizemode='area',
                sizeref=0.25*permax/100,
                showscale=True,
                colorbar={'title':"Mean",'len':1+1/len(selG)}),
                customdata=DOT['median'],
                hovertemplate='Group:%{x}<br>Gene:%{y}<br>Cell number:%{text}<br>Percentage:%{marker.size:.1f}<br>Mean:%{marker.color:.2f}<br>Median:%{customdata:.2f}'),
        row=subF['dot'][0],col=subF['dot'][1])
    perScale = [(int)(permax*i) for i in [1,0.67,0.33]]
    figDot.add_trace(go.Scatter(x=legendXY["x"],y=legendXY["y"],mode='markers+text',
            text=['%d%%'%i for i in perScale],
            textfont={'size':13,"color":"#000000"},
            marker=dict(size=[permax,permax*0.67,permax*0.33],
                color="#888888",
                sizemin=0,
                sizemode='area',
                sizeref=0.25*permax/100,
                showscale=False),
                hoverinfo='skip'),
        row=subF['legend'][0],col=subF['legend'][1])
    figDot.update_yaxes(linecolor="#000",showdividers=True,dividercolor="#444",title={"text":""},
                    tickfont={"size":15})
    figDot.update_xaxes(linecolor="#000",showdividers=True,dividercolor="#444",title={"text":""},
                    tickfont={"size":15},tickangle=45)
    figDot.update_xaxes(range=[-0.5,len(X[grp].cat.categories)],
                    row=subF['dot'][0],col=subF['dot'][1])
    figDot.update_xaxes(showticklabels=False,showline=False,fixedrange=True,
                        row=subF['legend'][0],col=subF['legend'][1])
    figDot.update_yaxes(showticklabels=False,showline=False,fixedrange=True,
                        row=subF['legend'][0],col=subF['legend'][1])
    if legendDirection=="h":
        figDot.update_annotations(x=0,y=1)#,xanchor="right"
        figDot.update_xaxes(range=[0.5,max(3.5,len(X[grp].cat.categories))],
                            row=subF['legend'][0],col=subF['legend'][1])
    else:
        figDot.update_yaxes(range=[min(0.5,5-len(selG)),3.5],
                            row=subF['legend'][0],col=subF['legend'][1])
    figDot.update_layout(plot_bgcolor='#fff',title={"text":grp},title_x=0.5,height=h,width=w,
                     margin={"l":0,"r":0,"t":30,"b":xlabH},showlegend=False)
    buffer = io.StringIO()
    figDot.write_html(buffer,include_plotlyjs=False,full_html=False)
    return buffer.getvalue()

def getAdditionalPara(argv):
    cN = None
    gCut = None
    expRange = [None,None]
    permax = None
    logmax = float('inf')
    lDirection = 'v'
    vname = None
    try:
        opts, args = getopt.getopt(argv,"n:g:l:e:p:d:v:",["ncell=","gcutoff=","logmax=","exprange=","percentagemax=","direction=","var="])
    except getopt.GetoptError:
        print("Usage: plotH5ad path/to/H5ad/file plot/type A/gene/list An/annotation/group -n cell/number -g gene/cutoff -l max/value/log -e min,max/exp/scale -p max/percentage/scale -d legend/direction -v <gene/column/header>")
        exit()
    for opt, arg in opts:
        if opt in ("-n", "--ncell"):
            cN = int(arg)
        elif opt in ("-g", "--gcutoff"):
            gCut = float(arg)
        elif opt in ("-l", "--logmax"):
            logmax = float(arg)
        elif opt in ("-e", "--exprange"):
            expRange = [float(x) for x in arg.split(",")]
            if len(expRange)!=2:
                print("expression scale format require min,max")
                exit()
        elif opt in ("-p", "--percentagemax"):
            permax = float(arg)
        elif opt in ("-d","--direction"):
            lDirection = arg[0].lower() # can be v, h, vertical, horizontal, V, H, V..., H...
            if lDirection not in ("v","h"):
                print("-d, --direction: only takes v, h, vertical, horizontal")
                exit()
        elif opt in ('-v',"-var"):
            vname = arg
    return cN,gCut,expRange,permax,logmax,lDirection,vname

def main():
    strH5ad = sys.argv[1]
    plotType = sys.argv[2]
    genes = sys.argv[3]
    grp = sys.argv[4]
    cN,gCut,expRange,permax,logmax,lDirection,vname = getAdditionalPara(sys.argv[5:])
    if plotType=='violin':
        print(violin(strH5ad,genes,grp,cN,gCut,vname))
    elif plotType=='dot':
        print(dot(strH5ad,genes,grp,cN,gCut,expRange,permax,logmax,lDirection,vname))
    else:
        print("ERROR: plot type ("+plotType+") is unknown!")
        exit()

if __name__ == "__main__":
    main()

#/share/cellxgene/VIP/bin/plotH5adPlotly.sh /raid/data/single_cell_portal/all_h5ad_files/f171db61-e57e-4535-a06a-35d8b6ef8f2b.h5ad violin A1BG 'cell_type' -n 10000 -g 0.1 -v feature_name
