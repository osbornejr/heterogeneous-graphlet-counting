using DataFrames

function tex_summary_table(dataframe::DataFrame,outfile::String)

    ## Take a raw count data table and make a summary table of it, including coding noncoding breakdown, how many transcripts 
    #separate breakdown of transcripts and samples?
    #TODO add check for :transcript_type column in dataframe.

    ##number of transcripts
    not = size(dataframe)[1]
    ##number of coding
    noc = size(filter(:transcript_type=>x->x=="coding",dataframe))[1]
    ##number of noncoding
    non = size(filter(:transcript_type=>x->x=="noncoding",dataframe))[1]

    ##number of samples
    

    table = "\\begin{tabular}"

    #number of rows and columns
    rows = 1 #size(dataframe)[1]
    columns = 5 #size(dataframe)[2]
    
    table *= "{@{}$(repeat("l",columns))@{}}\n"
    table *= "\\toprule\n"

    #add in column names
    for (i,n) in enumerate(["Number of transcripts","Coding transcripts","Coding percentage","Noncoding transcripts","Noncoding percentage"])
        #latex doesn't like _
        cn = replace(n,"_"=>"-")
        table *= "$(cn)"
        if (i<columns)
            table *=" & "
        end
    end
    table *= " \\\\\n"
    table *= "\\midrule\n"

    for (i,r) in enumerate([not,noc,(100*noc/not),non,(100*non/not)])
        for (j,d) in enumerate(r)
            #format floats and ints
            if (typeof(d)<:AbstractFloat)
                if (d>10000)
                    table *= " $(Int(round(d))) "
                else
                    table *= "$(round(d,digits=2)) "
                end
            else
                table *= "$(d) "
            end 
            if (j<columns)
                table *=" & "
            end
        end
        table *= "\\\\\n"
    end
    table*= "\\bottomrule\n"

    table*= "\\end{tabular}\n"

    write(outfile,table)
end
