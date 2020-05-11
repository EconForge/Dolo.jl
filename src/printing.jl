
sanitize(s, dynvars::Array{Symbol,1}) = s
sanitize(s::Symbol, dynvars::Array{Symbol,1})= s in dynvars ? :($s(0)) : s
function sanitize(eq::Expr, dynvars::Array{Symbol,1})
    if eq.head == :(=)
        return sanitize( :( $(eq.args[1])==$( eq.args[2] ) ), dynvars )
    else
        if eq.head == :call
            if eq.args[1] in dynvars
                return eq
            else
                return Expr(eq.head, [sanitize(e, dynvars) for e in eq.args]...)
            end
        else
            head = eq.head
            args = eq.args
            return Expr(eq.head, [sanitize(e, dynvars) for e in args]...)
        end
    end
end

function sanitize(a,model::Dolo.ASModel)
    dynvars = get_variables(model)
    return sanitize(a, dynvars)
end

function Base.show(io::IO, m::SModel{ID}) where ID
    println(io,
    """NumericModel
      - name: $(m.name)
      - filename: $(m.name)
    """)
end

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
