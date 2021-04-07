function Base.show(io::IO, ::MIME"text/html", model::Model)

    source = """
         <table>
             <td><b>Model</b></td>
         <tr>
            <td>name</td>
            <td>$(model.name)</td>
          </tr>
          <tr>
            <td>filename</td>
            <td>$(model.filename)</td>
          </tr>
        </table>
    """

    header = "<tr><td><b>Type</b></td><td><b>Equation</b></td></tr>\n"

    table_lines = []
    for eq_type in keys(model.equations)
        equations = model.equations[eq_type]
        if eq_type in [:controls_ub, :controls_lb]
            continue
        end
        for i=1:length(equations)
            if i==1
                eqt = eq_type
            else
                eqt = ""
            end
            eq = equations[i]
            eq = sanitize(eq, model)
            eqtex = Dolang.latex(eq)
            eqtex = replace(eqtex, "*"=>" ")
            fmt_line = "<tr><td>$eqt</td><td>\\[$eqtex\\]</td></tr>"
            push!(table_lines, fmt_line)
        end
    end

    source = string(source, "<table width=\"100%\">", header, table_lines..., "<table>")
    print(io, source)
end
