import Dolo

errs = Dolo.lint("LAMP_2s.yaml")

for err in errs
    println(err)
end
