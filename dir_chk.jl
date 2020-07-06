##
function folder_check(dirname)
    for i = 1:size(cfg.inputPlates,1)
        if(isdir("$(dirname)$(i)"))
            rm("$(dirname)$(i)", recursive=true)
        end
        mkdir("$(dirname)$(i)")
    end
end


function comln(fileId, var, val, comment = "1")
    if length(comment) > 1
        println(fileId, """#########################
        # $(comment)
        #########################""")
    end
    println(fileId, "$(var)\t$(val)")
end

function write_pdb(box, filename)
    fileId = open(filename, "w")
    println(fileId, "REMARK jectification generated")
    for i = 1:box.atomNum
        print(fileId, adj_left("ATOM",6))               #1-6
        print(fileId, adj_right(box.curAtom[i],5))      #7-11
        print(fileId, " ")                              #12whitespace
        print(fileId, adj_right(box.atomName[i],4))     #13-16
        print(fileId, " ")                              #17
        print(fileId, adj_right(box.molName[i],3))      #18-20
        print(fileId, " ")                              #21
        print(fileId, "A")                              #22
        print(fileId, adj_right(box.curMol[i],4))       #23-26
        print(fileId, " "^4)                            #27-30
        print(fileId, adj_right(box.x[i],8,3))          #31-38
        print(fileId, adj_right(box.y[i],8,3))          #39-46
        print(fileId, adj_right(box.z[i],8,3))          #47-54
        print(fileId, adj_right(1.0, 6, 2))             #55-60
        print(fileId, adj_right(0.0, 6, 2))             #61-66
        print(fileId," "^6)                             #67-72
        print(fileId, adj_left("M1", 4))                #73-76
        print(fileId, adj_right(box.element[i],2))      #77-78
        #print(fileId, adj_right(0.0,2))                 #79-80
        
        print(fileId, "\n")
    end
    print(fileId,"END")
    close(fileId)
end

function adj_left(var, s, dot=2)
    if typeof(var) == String
        var = lstrip(rstrip(var))
        if s > length(var)
            whadd = s - length(var)
            return var * " "^whadd
        else
            return var
        end
    elseif typeof(var) == Int
        des = floor(Int64, log10(var)) + 1
        if des < s
            return string(var) * " "^(s-des)
        else
            return string(var)
        end
    elseif typeof(var) == Float64
        return format(var, precision=s)
    end
end

function adj_right(var, s, dot = 2)
    if typeof(var) == String
        var = lstrip(rstrip(var))
        if s > length(var)
            whadd = s - length(var)
            return  " "^whadd * var
        else
            return var
        end
    elseif typeof(var) == Int
        des = floor(Int64, log10(var)) + 1
        if des < s
            return  " "^(s-des) * string(var)
        else
            return string(var)
        end
    elseif typeof(var) == Float64
        if var > 1 
            des = floor(Int64, log10(var)) + 1
        else
            des = 1
        end
        if des < s - dot-1
            return " "^(s-dot-des-1) * string(format(var, precision=dot))
        else
            return format(var, precision=dot)
        end
    end
end

