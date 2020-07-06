
function towhee_equlibrate(equlibrationCheck, equlibrationId)
#function towhee_equlibrate()    #debug editionion
    fileId = open("equlibration","r")
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) >1
        s = split(s,"+++++ end of markov chain +++++")[2]
    else
        println("towhee equlibration error")
        exit()
    end
    println(typeof(s),s)
    close(fileId)
end

a = 1
b = 1
towhee_equlibrate(a,b)
