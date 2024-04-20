module _WindowData
export WindowData
    """    
        WinData, Log = WindowData(Data) 

    Windows the sequence(s) given in `Data` into a collection of subsequnces  
    of floor(N/5) elements with no overlap, excluding any remainder
    elements that do not fill the final window.
    If `Data` is a univariate sequence (vector), `Windata` is a vector of 5
    vectors. If `Data` is a set of multivariate sequences (NxM matrix), 
    each of M columns is treated as a sequence with N elements 
    and `WinData` is a vector of 5 matrices of size [(floor*N,5), M]. 
    The `Log` dictionary contains information about the windowing process, including:
    `DataType`      - The type of data sequence passed as `Data`\n
    `DataLength`    - The number of sequence elements in `Data`\n
    `WindowLength`  - The number of elements in each window of `WinData`\n
    `WindowOverlap` - The number of overlapping elements between windows\n
    `TotalWindows`  - The number of windows extracted from `Data`\n
    `Mode`          - Decision to include or exclude any remaining sequence 
                        elements (`< WinLen`) that do not fill the window.\n

        WinData, Log = WindowData(Data::AbstractArray{T} where T<:Real, WinLen::Union{Nothing,Int}=nothing, Overlap::Int=0, Mode::String="exclude")

    Windows the sequence(s) given in `Data` into a collection of subsequnces
    using the specified keyword arguments: 

    # Arguments:
    `WinLen`  - Number of elements in each window, a positive integer (>10)\n
    `Overlap` - Number of overlapping elements between windows, a positive integer (< WinLen)\n
    `Mode`    - Decision to include or exclude any remaining sequence
                    elements (< `WinLen`) that do not fill the window,
                    a string - either `"include"` or `"exclude"` (default).\n

    # See also  `ExampleData`
    """
    function WindowData(Data::AbstractArray{T} where T<:Real;
                            WinLen::Union{Nothing,Int}=nothing, Overlap::Int=0, Mode::String="exclude")

    if ndims(Data)==1 
        DataType =  "single univariate vector (1 sequence)"
        N = length(Data)
        Dn = 0 
    elseif ndims(Data)==2 
        N, Dn = size(Data)
        DataType = "multivariate matrix ("*string(Dn)*" vectors)"
    else
        error("Only a vector or a Matrix can be passed as Data!")
    end
    (N>10) ? nothing : error("Data:   must be a numpy Vector (length N) or an NxM numpy matrix where N>10 and M>1")
    (isnothing(WinLen)) ? WinLen = Int(floor(N/5)) : nothing
    (10<WinLen<N) ? nothing : error("WinLen: must be an integer such that 10 < WinLen < N")
    (0<=Overlap<WinLen) ? nothing : error("Overlap: The number of overlapping window samples such that 0 <= Overlap < WinLen")
    (lowercase(Mode) in ["exclude","include"]) ? nothing : error("Mode: Option to include/exclude samples that do not fill final window, either 'exclude' or 'include'") 
            
    M = Int(floor((N - Overlap)/(WinLen - Overlap)))
    Step = Int(WinLen-Overlap)
    Xout = []
    map(k -> push!(Xout, Data[k:k+WinLen-1,:]), 1:Step:M*Step)        
    (lowercase(Mode)=="include" && (length(Xout)-1)*Step+WinLen!=N) ? (push!(Xout, Data[1+M*Step:end,:]); M+=1) : nothing

    Log = Dict("DataType" => DataType, "DataLength" => N, "WindowLength" => WinLen,
               "WindowOverlap" => Overlap,  "TotalWindows" => M,   "Mode" => Mode)
    
    return Xout, Log
    end
end

"""
    Copyright 2024 Matthew W. Flood, EntropyHub
    
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
    
     http://www.apache.org/licenses/LICENSE-2.0
    
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    
    For Terms of Use see https://github.com/MattWillFlood/EntropyHub
"""