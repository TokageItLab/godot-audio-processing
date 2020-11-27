class_name ComplexCalc

static func cprint(x: Vector2) -> String:
    var s: String = "" if x.x < 0 else "+"
    return String(x.x) + s + String(x.y) + "i"

static func cmlt(x1: Vector2, x2: Vector2) -> Vector2:
    var r: float = x1.x * x2.x - x1.y * x2.y
    var i: float = x1.x * x2.y + x1.y * x2.x
    return Vector2(r, i)

static func cdiv(x1: Vector2, x2: Vector2) -> Vector2:
    var r: float = x1.x * x2.x + x1.y * x2.y
    var i: float = x1.y * x2.x - x1.x * x2.y
    var d: float = x2.x * x2.x + x2.y * x2.y
    return Vector2(r / d, i / d)

static func cexp(x: Vector2) -> Vector2:
    var e = exp(x.x)
    return Vector2(e * cos(x.y), e * sin(x.y))

static func cpow(x: Vector2, n: float) -> Vector2:
    var r = sqrt(x.x * x.x + x.y * x.y)
    var t = atan(x.y / x.x)
    var p = pow(r, n)
    return Vector2(p * pow(cos(t), n), p * pow(sin(t), n))

static func cabs(x: Vector2) -> float:
    return x.length()
