
function RTs_SS(a,b,distri,cohmeasure)
    # a quantity to mean
    # b quantity to test
    # Return the mean of RTs and STD
    # Sucess after Sucess
    output = Float64[]
    for i=1:length(a)-1
        x = b[i]
        y=b[i+1]
        x > 0 && y >0 && distri[i+1]*distri[i] > 0 && distri[i+1] == cohmeasure || continue
        push!(output,a[i+1])
    end
    return mean(output),std(output),length(output)
end

function performance(results)
    # compute the performance of the model
    # need the list of of results for a specific coherence
    return sum(results)/length(results)
end

function RTs_SE(a,b,distri,cohmeasure)
    # a quantity to mean
    # b quantity to test
    # Error after Sucess
    t = 0
    c = 0
    output=Float64[]
    for i=1:length(a)-1
        x = b[i]
        y=b[i+1]
        x > 0 && y <0 && distri[i+1]*distri[i] > 0 && distri[i+1] == cohmeasure || continue
        push!(output,a[i+1])

    end
    return mean(output),std(output),length(output)
end

function RTs_ES(a,b,distri,cohmeasure)
    # a quantity to mean
    # b quantity to test
    ## Sucess after Error
    output=Float64[]
    for i=1:length(a)-1
        x = b[i]
        y=b[i+1]
        x < 0 && y >0 && distri[i+1] == cohmeasure|| continue

        push!(output,a[i+1])

    end
    return mean(output),std(output),length(output)
end

function RTs_SSS(a,b)
    # Sucess after Sucess after Sucess
    # same sign coh but not same one
    output=Float64[]
    for i=1:length(a)-2
        x = b[i]
        y=b[i+1]
        z=b[i+2]
        x > 0 && y >0 && z>0 && distri[i]*distri[i+1]>0 && distri[i+1]*distri[i+2] && distri[i+2] == cohmeasure|| continue
        push!(output,a[i+2])
    end
    return mean(output),std(output),length(output)
end


function RTs_SSE(a,b)
    # Error after Sucess after Sucess
    t = 0
    c = 0
    output=Float64[]

    for i=1:length(a)-2
        x = b[i]
        y=b[i+1]
        z=b[i+2]
        x > 0 && y >0 && z<0&& distri[i]*distri[i+1]>0 && distri[i+1]*distri[i+2] && distri[i+2] == cohmeasure|| continue
        push!(output,a[i+2])
    end
    return mean(output),std(output),length(output)
end

function RTs_SS_chgt(RT,Res,distri,coherence)
    output = Float64[]
    for i =1:length(RT)-1
        if Res[i]==1 && Res[i+1] == 1 && distri[i+1] == coherence && distri[i]*distri[i+1] <0  || continue
            push!(output,RT[i+1])
        end
    end
    return mean(output)
end

function analysis_RTs_E(result,Rts_coh)
        listc=Float64[]
        listv=Float64[]

        for c in [0 1 2 3 4 5 6 7 8]
                tempv = Float64[]
                temp=Float64[]
                for i in 1:length(result["distri"])-1
                        if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                                result["Confidence"][i+1]>(c-1) &&result["NotError"][i]==-1 #&&result["NotError"][i+1]==1
                                x=result["distri"][i+1]
                                push!(tempv,result["RTs"][i+1])
                                if abs(x) <10
                                push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
end
                        end
                end
                push!(listc,mean(temp))
                push!(listv,var(temp))
        end
        return listc,listv
end

function analysis_RTs_SE(result,Rts_coh)
    listc=Float64[]
    listv=Float64[]

    for c in [0 1 2 3 4 5 6 7 8]
            tempv = Float64[]
            temp=Float64[]
            for i in 1:length(result["distri"])-1
                    if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                            result["Confidence"][i+1]>(c-1) &&result["NotError"][i]==-1 &&result["NotError"][i+1]==1
                            x=result["distri"][i+1]
                            push!(tempv,result["RTs"][i+1])
                            push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
                    end
            end
            push!(listc,mean(temp))
            push!(listv,var(temp))
    end
    return listc,listv
end

function analysis_RTs_SS(result,Rts_coh)
    listc=Float64[]
    listv=Float64[]

    for c in [0 1 2 3 4 5 6 7 8]
            tempv = Float64[]
            temp=Float64[]
            for i in 1:length(result["distri"])-1
                    if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                            result["Confidence"][i+1]>(c-1) &&result["NotError"][i]==1 &&result["NotError"][i+1]==1
                            x=result["distri"][i+1]
                            push!(tempv,result["RTs"][i+1])
                            push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
                    end
            end
            push!(listc,mean(temp))
            push!(listv,var(temp))
    end
    return listc,listv
end

function analysis_RTs_E_cohspec(result,Rts_coh)
        listc=Float64[]
        for c in [0 1 2 3 4 5 6 7 8]
                temp=Float64[]
                for i in 1:length(result["distri"])-1
                        if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                                result["Confidence"][i+1]==c &&result["NotError"][i]==-1 #&&result["NotError"][i+1]==1
                                x=result["distri"][i+1]
                                push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
                        end
                end
                push!(listc,mean(temp))
        end
        return listc
end


function analysis_RTs_SE_cohspec(result,Rts_coh)
        listc=Float64[]
        for c in [0 1 2 3 4 5 6 7 8]
                temp=Float64[]
                for i in 1:length(result["distri"])-1
                        if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                                result["Confidence"][i+1]==c &&result["NotError"][i]==-1 &&result["NotError"][i+1]==1
                                x=result["distri"][i+1]
                                push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
                        end
                end
                push!(listc,mean(temp))
        end
        return listc
end


function analysis_RTs_SS_cohspec(result,Rts_coh)
        listc=Float64[]
        for c in [0 1 2 3 4 5 6 7 8]
                temp=Float64[]
                for i in 1:length(result["distri"])-1
                        if sign(result["distri"][i])==sign(result["distri"][i+1]) &&
                                result["Confidence"][i+1]==(c-1) &&result["NotError"][i]==1 &&result["NotError"][i+1]==1
                                x=result["distri"][i+1]
                                push!(temp,(Rts_coh["$x"]-result["RTs"][i+1])/Rts_coh["$x"])
                        end
                end
                push!(listc,mean(temp))
        end
        return listc
end

function load_Jerome_E()

       data = CSV.read("/home/kevin/Documents/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
        #data = CSV.read("/media/kevin/Boulot/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)

        rt_confidence = [] # for confidence task
        results_confidence = []
        cursor_confidence = []
        coherence_confidence = []
        reponse_confidence = []
        name_confidence=[]
        number=Float64[]

        for i=1:length(data[1])
                if  data[3][i] == "confidence"
                        push!(reponse_confidence,data[6][i])
                        push!(cursor_confidence,data[11][i])
                        push!(coherence_confidence,data[5][i])
                        push!(rt_confidence,data[7][i])
                        push!(results_confidence,data[20][i])
                        push!(name_confidence,data[2][i])
                        push!(number,parse(Float64,data[1][i]))
                    end
        end
        #list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        #list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4"]
        list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4" ]

        dic_rt_confidence = Dict()
        for coh_string in list_coherence_string
                rt_new = rt_confidence[coherence_confidence.==parse(Float64,coh_string)]
                dic_rt_confidence[coh_string]=mean(rt_new[rt_new.<1500])
        end
        dic_rt_cursor_error=Dict()
        for cursor in [0 1 2 3 4 5 6 7 8]
             rt_list_difference=Float64[]

             for i=1:(length(rt_confidence)-1)
                if  sign(coherence_confidence[i])==sign(coherence_confidence[i+1])  &&
                    name_confidence[i] == name_confidence[i+1] &&
                    results_confidence[i] == 0 && number[i] +1 == number[i+1]#&&results_confidence[i+1] == 1

                   if rt_confidence[i]<1500 && rt_confidence[i+1] <1500 &&
                      cursor_confidence[i]>cursor #&& cohVEA[i]!=0


                        x=coherence_confidence[i+1]
                        if x==0 # because 0 is a float in our data
                            push!(rt_list_difference,((dic_rt_confidence["0"]-rt_confidence[i+1])/dic_rt_confidence["0"]))
                        else
                            push!(rt_list_difference,((dic_rt_confidence["$x"]-rt_confidence[i+1])/dic_rt_confidence["$x"]))

                    end
                    end
                end
             end
    println(length(rt_list_difference))
    dic_rt_cursor_error["$cursor"]=mean(rt_list_difference)*100
    dic_rt_cursor_error["var $cursor"]=var(rt_list_difference*100)/700
    dic_rt_cursor_error["label $cursor"]=length(rt_list_difference)
    end

    exp_list2=Float64[]
    var_list2=Float64[]
    for cursor in [0 1 2 3 4 5 6 7 8]
                push!(exp_list2,dic_rt_cursor_error["$cursor"])
                push!(var_list2,dic_rt_cursor_error["var $cursor"])
    end
    return exp_list2,var_list2
end

function load_Jerome_SS()

            data = CSV.read("/home/kevin/Documents/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
        #    data = CSV.read("/media/kevin/Boulot/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
            rt_confidence = [] # for confidence task
        results_confidence = []
        cursor_confidence = []
        coherence_confidence = []
        reponse_confidence = []
        name_confidence=[]
        number=[]
        for i=1:length(data[1])
                if  data[3][i] == "confidence"
                        push!(reponse_confidence,data[6][i])
                        push!(cursor_confidence,data[11][i])
                        push!(coherence_confidence,data[5][i])
                        push!(rt_confidence,data[7][i])
                        push!(results_confidence,data[20][i])
                        push!(name_confidence,data[2][i])
                        push!(number,parse(Float64,data[1][i]))
                    end
        end
        list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4"]
        dic_rt_confidence = Dict()
        for coh_string in list_coherence_string
                rt_new = rt_confidence[coherence_confidence.==parse(Float64,coh_string)]
                dic_rt_confidence[coh_string]=mean(rt_new[rt_new.<1500])
                println(mean(rt_new[rt_new.<1500]))
        end
        dic_rt_cursor_error=Dict()
        for cursor in [0 1 2 3 4 5 6 7 8]
             rt_list_difference=Float64[]

             for i=1:(length(rt_confidence)-1)
                if  coherence_confidence[i]==coherence_confidence[i+1]  &&
                    name_confidence[i] == name_confidence[i+1] &&
                    number[i]+1==number[i+1] &&
                    results_confidence[i] == 1 &&results_confidence[i+1] == 1

                   if rt_confidence[i]<1500 && rt_confidence[i+1] <1500 &&
                      cursor_confidence[i]>cursor #&& cohVEA[i]!=0


                        x=coherence_confidence[i+1]
                        if x==0 # because 0 is a float in our data
                            push!(rt_list_difference,((rt_confidence[i]-rt_confidence[i+1])/dic_rt_confidence["0"]))
                        else
                            push!(rt_list_difference,((rt_confidence[i]-rt_confidence[i+1])/dic_rt_confidence["$x"]))
                        end

                        #if x==0 # because 0 is a float in our data
                        #    push!(rt_list_difference,((dic_rt_confidence["0"]-rt_confidence[i+1])/dic_rt_confidence["0"]))
                        #else
                        #    push!(rt_list_difference,((dic_rt_confidence["$x"]-rt_confidence[i+1])/dic_rt_confidence["$x"]))

                    #end
                    end
                end
             end
    println(length(rt_list_difference))
    dic_rt_cursor_error["$cursor"]=mean(rt_list_difference)*100
    dic_rt_cursor_error["var $cursor"]=var(rt_list_difference)
    dic_rt_cursor_error["label $cursor"]=length(rt_list_difference)
    end

    exp_list2=Float64[]
    var_list2=Float64[]
    for cursor in [0 1 2 3 4 5 6 7 8]
                push!(exp_list2,dic_rt_cursor_error["$cursor"])
                push!(var_list2,dic_rt_cursor_error["var $cursor"])
    end
    return exp_list2,var_list2
end


function load_Jerome_SE_cohspec()

            data = CSV.read("/home/kevin/Documents/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
        #    data = CSV.read("/media/kevin/Boulot/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
    rt_confidence = [] # for confidence task
        results_confidence = []
        cursor_confidence = []
        coherence_confidence = []
        reponse_confidence = []
        name_confidence=[]

        for i=1:length(data[1])
                if  data[3][i] == "confidence"
                        push!(reponse_confidence,data[6][i])
                        push!(cursor_confidence,data[11][i])
                        push!(coherence_confidence,data[5][i])
                        push!(rt_confidence,data[7][i])
                        push!(results_confidence,data[20][i])
                        push!(name_confidence,data[2][i])
                    end
        end
        list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4"]
        dic_rt_confidence = Dict()
        for coh_string in list_coherence_string
                rt_new = rt_confidence[coherence_confidence.==parse(Float64,coh_string)]
                dic_rt_confidence[coh_string]=mean(rt_new[rt_new.<1500])
        end
        dic_rt_cursor_error=Dict()
        for cursor in [0 1 2 3 4 5 6 7 8]
             rt_list_difference=Float64[]

             for i=1:(length(rt_confidence)-1)
                if  sign(coherence_confidence[i])==sign(coherence_confidence[i+1])  &&
                    name_confidence[i] == name_confidence[i+1] &&
                    results_confidence[i] == 0 &&results_confidence[i+1] == 1

                   if rt_confidence[i]<1500 && rt_confidence[i+1] <1500 &&
                      cursor_confidence[i]==cursor+1 #&& cohVEA[i]!=0


                        x=coherence_confidence[i+1]
                        if x==0 # because 0 is a float in our data
                            push!(rt_list_difference,((dic_rt_confidence["0"]-rt_confidence[i+1])/dic_rt_confidence["0"]))
                        else
                            push!(rt_list_difference,((dic_rt_confidence["$x"]-rt_confidence[i+1])/dic_rt_confidence["$x"]))
                        end
                    end
                end
             end
    println(length(rt_list_difference))
    dic_rt_cursor_error["$cursor"]=mean(rt_list_difference)*100
    dic_rt_cursor_error["var $cursor"]=var(rt_list_difference)
    dic_rt_cursor_error["label $cursor"]=length(rt_list_difference)
    end

    exp_list2=Float64[]
    var_list2=Float64[]
    for cursor in [0 1 2 3 4 5 6 7 8]
                push!(exp_list2,dic_rt_cursor_error["$cursor"])
                push!(var_list2,dic_rt_cursor_error["var $cursor"])
    end
    return exp_list2,var_list2
end

function load_Jerome_SE()

            data = CSV.read("/home/kevin/Documents/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
        #    data = CSV.read("/media/kevin/Boulot/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
rt_confidence = [] # for confidence task
        results_confidence = []
        cursor_confidence = []
        coherence_confidence = []
        reponse_confidence = []
        name_confidence=[]
        number=Float64[]
        for i=1:length(data[1])
                if  data[3][i] == "confidence"
                        push!(reponse_confidence,data[6][i])
                        push!(cursor_confidence,data[11][i])
                        push!(coherence_confidence,data[5][i])
                        push!(rt_confidence,data[7][i])
                        push!(results_confidence,data[20][i])
                        push!(name_confidence,data[2][i])
                        push!(number,parse(Float64,data[1][i]))
                    end
        end
        list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4"]
        dic_rt_confidence = Dict()
        for coh_string in list_coherence_string
                rt_new = rt_confidence[coherence_confidence.==parse(Float64,coh_string)]
                dic_rt_confidence[coh_string]=mean(rt_new[rt_new.<1500])
        end
        dic_rt_cursor_error=Dict()
        for cursor in [0 1 2 3 4 5 6 7 8]
             rt_list_difference=Float64[]

             for i=1:(length(rt_confidence)-1)
                if  sign(coherence_confidence[i])==sign(coherence_confidence[i+1])  &&
                    name_confidence[i] == name_confidence[i+1] &&
                    number[i] + 1 == number[i+1] &&
                    results_confidence[i] == 0 &&results_confidence[i+1] == 1

                   if rt_confidence[i]<1500 && rt_confidence[i+1] <1500 &&
                      cursor_confidence[i]>=cursor+1 #&& cohVEA[i]!=0


                        x=coherence_confidence[i+1]
                        if x==0 # because 0 is a float in our data
                            push!(rt_list_difference,((dic_rt_confidence["0"]-rt_confidence[i+1])/dic_rt_confidence["0"]))
                        else
                            push!(rt_list_difference,((dic_rt_confidence["$x"]-rt_confidence[i+1])/dic_rt_confidence["$x"]))
                        end
                    end
                end
             end
    println(length(rt_list_difference))
    dic_rt_cursor_error["$cursor"]=mean(rt_list_difference)*100
    dic_rt_cursor_error["var $cursor"]=var(rt_list_difference)
    dic_rt_cursor_error["label $cursor"]=length(rt_list_difference)
    end

    exp_list2=Float64[]
    var_list2=Float64[]
    for cursor in [0 1 2 3 4 5 6 7 8]
                push!(exp_list2,dic_rt_cursor_error["$cursor"])
                push!(var_list2,dic_rt_cursor_error["var $cursor"])
    end
    return exp_list2,var_list2
end


function load_Jerome_nextRTE()

            data = CSV.read("/home/kevin/Documents/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
        #    data = CSV.read("/media/kevin/Boulot/PhD/PhD_Data/Jerome_17/secondFirst.csv",nullable=false)
rt_confidence = [] # for confidence task
        results_confidence = []
        cursor_confidence = []
        coherence_confidence = []
        reponse_confidence = []
        name_confidence=[]
        number=Float64[]
        for i=1:length(data[1])
                if  data[3][i] == "confidence"
                        push!(reponse_confidence,data[6][i])
                        push!(cursor_confidence,data[11][i])
                        push!(coherence_confidence,data[5][i])
                        push!(rt_confidence,data[7][i])
                        push!(results_confidence,data[20][i])
                        push!(name_confidence,data[2][i])
                        push!(number,parse(Float64,data[1][i]))
                    end
        end
        list_coherence = [-2.4 -1.6 -0.8 0 0.8 1.6 2.4]
        list_coherence_string=["-2.4" "-1.6" "-0.8" "0" "0.8" "1.6" "2.4"]
        dic_rt_confidence = Dict()
        for coh_string in list_coherence_string
                rt_new = rt_confidence[coherence_confidence.==parse(Float64,coh_string)]
                dic_rt_confidence[coh_string]=mean(rt_new[rt_new.<1500])
        end
        dic_rt_cursor_error=Dict()
        for cursor in [0 1 2 3 4 5 6 7 8]
             rt_list_difference=Float64[]

             for i=1:(length(rt_confidence)-1)
                if  name_confidence[i] == name_confidence[i+1] &&
                    number[i] + 1 == number[i+1] &&
                    results_confidence[i] == 0 # &&results_confidence[i+1] == 1

                   if rt_confidence[i]<1500 && rt_confidence[i+1] <1500 &&
                      cursor_confidence[i]==cursor+1 #&& cohVEA[i]!=0


                        x=coherence_confidence[i+1]
                        if x==0 # because 0 is a float in our data
                            push!(rt_list_difference,((dic_rt_confidence["0"]-rt_confidence[i+1])))
                        else
                            push!(rt_list_difference,((dic_rt_confidence["$x"]-rt_confidence[i+1])))
                        end
                    end
                end
             end
    println(length(rt_list_difference))
    dic_rt_cursor_error["$cursor"]=mean(rt_list_difference)
    dic_rt_cursor_error["var $cursor"]=var(rt_list_difference)
    dic_rt_cursor_error["label $cursor"]=length(rt_list_difference)
    end

    exp_list2=Float64[]
    var_list2=Float64[]
    for cursor in [0 1 2 3 4 5 6 7 8]
                push!(exp_list2,dic_rt_cursor_error["$cursor"])
                push!(var_list2,dic_rt_cursor_error["var $cursor"])
    end
    return exp_list2,var_list2
end
