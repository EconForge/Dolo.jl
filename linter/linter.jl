__precompile__(true)

module DoloLinter

    import YAML

    function yaml_node_from_file(fn::AbstractString)

        # Didn't find how to access the top node more easily.
        #
        yml_types = Dict{AbstractString,Function}()
        yml_types["!Cartesian"] = ((c,n)->n)
        yml_types["!Smolyak"] = ((c,n)->n)
        yml_types["!Normal"] = ((c,n)->n)
        yml_types["!MarkovChain"] = ((c,n)->n)
        yml_types["!Product"] = ((c,n)->n)
        yml_types["!PoissonProcess"] = ((c,n)->n)
        yml_types["!DeathProcess"] = ((c,n)->n)
        yml_types["!AgingProcess"] = ((c,n)->n)
        yml_types["!VAR1"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:null"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:bool"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:int"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:float"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:binary"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:timestamp"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:omap"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:pairs"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:set"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:str"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:seq"] = ((c,n)->n)
        yml_types["tag:yaml.org,2002:map"] = ((c,n)->n)


        return YAML.load_file(fn, yml_types)
    end

    Base.getindex(d::YAML.SequenceNode, k::Integer) = d.value[k]
    Base.keys(s::YAML.MappingNode) =  [e[1].value for e=s.value]
    Base.values(s::YAML.MappingNode) =  [e[2].value for e=s.value]
    Base.getindex(d::YAML.MappingNode, s::AbstractString) = d.value[findfirst(keys(d),s)][2]
    Base.getindex(d::YAML.MappingNode, s::Symbol) = d[string(s)]
    # extended index: find


    type Location
        start_mark::YAML.Mark
        end_mark::YAML.Mark
    end

    type MissingElement <: Exception
        path::Vector{AbstractString}
        k::Integer # index of first missing element
        loc::Location
    end

    function Base.getindex(d::YAML.MappingNode, l::AbstractString...)
        el = d
        for k=1:length(l)
            # if !(typeof(el)<:YAML.MappingNode)
            #     println("This is not going to fly")
            if !(l[k] in keys(el))
                loc = Location(el.start_mark, el.end_mark)
                throw(MissingElement(collect(l), k, loc))
            else
                el = el[l[k]]
            end
        end
        return el
    end
    function Base.getindex(d::YAML.MappingNode, l::Symbol...)
        return Base.getindex(d, [string(e) for e in l]...)
    end


    type LinterException <: Exception
        msg::AbstractString
        loc::Location
        src::AbstractString
    end

    type LinterWarning <: Exception
        msg::AbstractString
        loc::Location
        src::AbstractString
    end

    LinterException(msg::AbstractString, loc::Location) = LinterException(msg, loc, "<string>")

    LinterWarning(msg::AbstractString, loc::Location) = LinterWarning(msg, loc, "<string>")

    function get_loc(d)
        Location(d.start_mark, d.end_mark)
    end

    function check(d::YAML.MappingNode, filename="<string>")

        ### so far, this is just an example

        errors = []
        warnings = []

        # check symbol names:
        required_symbol_types = ["states", "controls", "exogenous", "parameters"]
        optional_symbol_types = [ "values", "rewards", "expectations"]
        known_symbol_types = cat(1, required_symbol_types, optional_symbol_types)
        sym_names = keys(d[:symbols])
        for s in required_symbol_types
            if !(s in sym_names)
                msg = string("Missing symbol of type: ", s)
                push!(errors, LinterWarning(msg, get_loc(d[:symbols]), filename))
            end
        end
        for (i,sg) in enumerate(keys(d["symbols"]))
            if !(sg in cat(1,known_symbol_types))
                msg = string("Unknown symbol type: ", sg)
                # get precise location
                loc = get_loc(d["symbols"].value[i][1])
                push!(warnings, LinterWarning(msg, loc, filename))
            end
        end

        return [errors, warnings]

    end

    function check(fn::AbstractString)
        d = yaml_node_from_file(fn)
        errs, wars = check(d, fn)
    end

end
