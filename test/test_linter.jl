import Dolo

# an example model shouldn't raise any error
models_dir = joinpath(Dolo.pkg_path,"examples","models")
for fname in (readdir(models_dir))
    filename = joinpath(models_dir, fname)
    if isfile(filename)
        println(filename)
        ch = Dolo.lint(filename)
        @assert length(ch) == 0
    end
end

# for now, we check that problematic models don't make linter failxx
errors_dir = joinpath(Dolo.pkg_path,"linter","models")
for fname in (readdir(errors_dir))
    filename = joinpath(errors_dir, fname)
    if isfile(filename)
        ok = true
        try
            ch = Dolo.lint(filename)
        catch err
            ok = false
            println("Uncaught error while parsing")
            println(err)
        end
        @assert ok
    end
end
