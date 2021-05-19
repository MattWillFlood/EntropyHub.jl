module _ExampleData
export ExampleData
using DelimitedFiles
using HTTP

    function ExampleData(SigName)

    Chk = ["uniform","uniform2","randintegers2","randintegers",
        "henon","chirp","gaussian","gaussian2","lorenz",
        "uniform_Mat", "gaussian_Mat", "entropyhub_Mat",
        "mandelbrot_Mat","randintegers_Mat"]

    SigName in Chk ? nothing : error("SigName must be one of the following:\n$Chk")
    url = "https://raw.githubusercontent.com/MattWillFlood/EntropyHub/main/ExampleData/" * SigName * ".txt"
    Temp = HTTP.get(url).body

        if (SigName in ["henon","lorenz"]) || (SigName[end] == '2') || (SigName[end-2:end] == "Mat")  
            X = readdlm(Temp, skipstart=2)    
        else
            X = readdlm(Temp, skipstart=2);
            size(X)[1] < size(X)[2] ? X = X[1,:] : X = X[:,1]
        end

    return X
    end

end