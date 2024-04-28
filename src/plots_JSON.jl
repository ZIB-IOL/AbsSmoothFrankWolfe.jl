using JSON

function save_to_json(example, iterates_information_vector)
    dict = Dict{Symbol, Any}()
    dict[:iteration] = [i[1] for i in iterates_information_vector]
    dict[:x] = [i[2] for i in iterates_information_vector]
    dict[:objective_value] = [i[3] for i in iterates_information_vector]

    iter = dict[:iteration][end]
    file_name = "json/" * example * "_" * string(iter) * ".json"
    open(file_name, "w") do f
        JSON.print(f, dict) 
    end

end

function save_to_json2(example, iterates_information_vector)
    dict = Dict{Symbol, Any}()
    dict[:iteration] = [i[1] for i in iterates_information_vector]
    dict[:x] = [i[2] for i in iterates_information_vector]
    dict[:objective_value] = [i[3] for i in iterates_information_vector]
    dict[:dual_gap] = [i[4] for i in iterates_information_vector]

    iter = dict[:iteration][end]
    file_name = "json/" * example * "_" * string(iter) * ".json"
    open(file_name, "w") do f
        JSON.print(f, dict) 
    end

end
