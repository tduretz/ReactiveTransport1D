using Unitful

°C = u"°C"
°F = u"°F"
μm = 1e-6u"m"
m  = u"m"
Ra = u"Ra"
K  = u"K"
F  = u"kg"*u"m"^2/u"s"
mm = 1e-3u"m"

@show uconvert(°C, 212°F)

@show upreferred(F/m)

