import glob

#noCuts_njets = glob.glob("../outputPlots/h_njets.pdf")
#noCuts_ntops = glob.glob("../outputPlots/h_ntops.pdf")
#noCuts_nbjet = glob.glob("../outputPlots/h_nb.pdf")
#noCuts_HT = glob.glob("../outputPlots/h_ht.pdf")
#noCuts_HT = glob.glob("../outputPlots/h_fisher.pdf")
#noCuts_HT = glob.glob("../outputPlots/h_bdt.pdf")
noCuts_njets = glob.glob("../outputPlots/h_njets.pdf")
noCuts_ntops = glob.glob("../outputPlots/blind_ntops.pdf")
noCuts_nbjet = glob.glob("../outputPlots/blind_nb.pdf")
noCuts_HT = glob.glob("../outputPlots/blind_ht.pdf")
noCuts_HT = glob.glob("../outputPlots/blind_fisher.pdf")
noCuts_HT = glob.glob("../outputPlots/blind_bdt.pdf")

#njets_stack_plots_0l = glob.glob("../outputPlots/h_njets_0l_ge6j_HT500_ge2b*.pdf")
#ntops_stack_plots_0l = glob.glob("../outputPlots/h_ntops_0l_ge6j_HT500_ge2b*.pdf")
#nb_stack_plots_0l = glob.glob("../outputPlots/h_nb_0l_ge6j_HT500_ge2b*.pdf")
#ht_stack_plots_0l = glob.glob("../outputPlots/h_HT_0l_ge6j_HT500_ge2b*.pdf")
#fisher_stack_plots_0l = glob.glob("../outputPlots/h_fisher_0l_ge6j_HT500_ge2b*.pdf")
njets_stack_plots_0l = glob.glob("../outputPlots/h_njets_0l_ge6j_HT500_ge2b*.pdf")
ntops_stack_plots_0l = glob.glob("../outputPlots/blind_ntops_0l_ge6j_HT500_ge2b*.pdf")
nb_stack_plots_0l = glob.glob("../outputPlots/blind_nb_0l_ge6j_HT500_ge2b*.pdf")
ht_stack_plots_0l = glob.glob("../outputPlots/blind_HT_0l_ge6j_HT500_ge2b*.pdf")
fisher_stack_plots_0l = glob.glob("../outputPlots/blind_fisher_0l_ge6j_HT500_ge2b*.pdf")

fisher_plots_0l = glob.glob("../outputPlots/fisher_*.pdf") 
fisherNorm_plots_0l = glob.glob("../outputPlots/fisherNorm_*.pdf") 
fisherRocCompare_plots_0l = glob.glob("../outputPlots/fisherRocCompare_*.pdf") 
fisherRoc_plots_0l = glob.glob("../outputPlots/fisherRoc_*.pdf") 

plotList = [noCuts_njets,noCuts_ntops,noCuts_nbjet,noCuts_HT]

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

def write(f,glob,title):
    for i in range(0, len(glob), 2):
        if i == len(glob) - 1:
            f.write(tex_snippet_1plot % (title, glob[i].replace(".pdf","")))
            print glob[i]
        else:
            f.write(tex_snippet % (title, glob[i].replace(".pdf",""), glob[i+1].replace(".pdf","")))
            print glob[i]
            print glob[i+1]

def writeList(f,l,title):
    for i in range(0, len(l), 2):
        if i == len(l) - 1:
            f.write(tex_snippet_1plot % (title, l[i][0].replace(".pdf","")))
            print l[i][0]
        else:
            f.write(tex_snippet % (title, l[i][0].replace(".pdf",""), l[i+1][0].replace(".pdf","")))
            print l[i][0]
            print l[i+1][0]

            
f = open("stack_snippet.tex",'w')    
write(f,fisher_plots_0l,"0 lepton - Fisher plots")
#writeList(f,plotList,"No cuts")
write(f,fisherNorm_plots_0l,"0 lepton - Normalized Fisher plots")
write(f,fisherRocCompare_plots_0l,"0 lepton - Fisher Roc plots: Test vs V3")
write(f,fisherRoc_plots_0l,"0 lepton - Fisher Roc plots: Test")
write(f,njets_stack_plots_0l,"0 lepton - stack plots")
write(f,ntops_stack_plots_0l,"0 lepton - stack plots")
write(f,nb_stack_plots_0l,"0 lepton - stack plots")
write(f,ht_stack_plots_0l,"0 lepton - stack plots")
write(f,fisher_stack_plots_0l,"0 lepton - stack plots")

f.close()
