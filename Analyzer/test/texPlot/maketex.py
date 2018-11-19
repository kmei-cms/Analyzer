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
    path = "../outputPlots/"
    f = open("stack_snippet.tex",'w')

    write(f, path+"blind_njets*.pdf",   "NJet - stack plots")
    write(f, path+"blind_ntops*.pdf",   "NTops - stack plots")
    write(f, path+"blind_nb*.pdf",      "NBottoms - stack plots")
    write(f, path+"blind_deepESM*.pdf", "DeepESM - stack plots")
    write(f, path+"blind_ht*.pdf",      "HT - stack plots")
    write(f, path+"blind_lPt*.pdf",     "Lepton PT - stack plots")
    write(f, path+"blind_lEta*.pdf",    "Lepton Eta - stack plots")
    write(f, path+"blind_mbl*.pdf",     "mbl - stack plots")

    f.close()

if __name__ == '__main__':
    main()
