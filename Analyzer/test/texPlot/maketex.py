import glob

tex_snippet = """%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    \\begin{frame}{%s}
    
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
    
    \\end{frame}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """
    
tex_snippet_1plot = """%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    \\begin{frame}{%s}
    
        \\includegraphics[width=0.5\\textwidth]{{"%s"}.pdf}
    
    \\end{frame}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """

def write(f,globString,title):
    plotList = glob.glob(globString)

    for i in range(0, len(plotList), 2):
        if i == len(plotList) - 1:
            f.write(tex_snippet_1plot % (title, plotList[i].replace(".pdf","")))
            print plotList[i]
        else:
            f.write(tex_snippet % (title, plotList[i].replace(".pdf",""), plotList[i+1].replace(".pdf","")))
            print plotList[i]
            print plotList[i+1]

def main():    
    paths = [("../outputPlots/FullRun2_2017/","2017"), ("../outputPlots/FullRun2_2018/","2018")]
    f = open("stack_snippet.tex",'w')

    for path in paths:
        write(f, path[0]+"h_njets*.pdf",       path[1]+" NJet - stack plots")
        write(f, path[0]+"blind_ntops*.pdf",   path[1]+" NTops - stack plots")
        write(f, path[0]+"blind_nb*.pdf",      path[1]+" NBottoms - stack plots")
        write(f, path[0]+"blind_deepESM*.pdf", path[1]+" DeepESM - stack plots")
        write(f, path[0]+"blind_deepESMMerged*.pdf", path[1]+" DeepESMMerged - stack plots")
        write(f, path[0]+"blind_ht*.pdf",      path[1]+" HT - stack plots")
        write(f, path[0]+"blind_lPt*.pdf",     path[1]+" Lepton PT - stack plots")
        write(f, path[0]+"blind_lEta*.pdf",    path[1]+" Lepton Eta - stack plots")
        write(f, path[0]+"blind_mbl*.pdf",     path[1]+" mbl - stack plots")
        write(f, path[0]+"fisherNorm_h_deepESM*.pdf", path[1]+" Norm DeepESM")
        write(f, path[0]+"fisherNorm_h_njets*.pdf",   path[1]+" Norm NJets")
        write(f, path[0]+"fisherRocCompare*.pdf",     path[1]+" Roc Curves")

    f.close()

if __name__ == '__main__':
    main()
