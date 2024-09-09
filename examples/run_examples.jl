
specific_files = [
    "small_example.jl",  
    "large_example.jl"
]

current = abspath(@__FILE__)


for file in specific_files
    full_path = abspath(file)
    
    if full_path == current
        println("Skipping the example: $file")
        continue
    end
    
    if !isfile(full_path)
        println("File not found: $file")
        continue
    end

    println("Running: $file")
    try
        include(full_path)  
    catch e
        println("Error running $file: $e")
    end
end


GC.gc()

