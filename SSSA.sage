p = 16857450949524777441941817393974784044780411511252189319
A = 16857450949524777441941817393974784044780411507861094535
B = 77986137112576

# y^2 = x^3 + A * x + B
tE = EllipticCurve(GF(p), [A, B])

tP = tE(5732560139258194764535999929325388041568732716579308775,
        14532336890195013837874850588152996214121327870156054248)
tQ = tE(2609506039090139098835068603396546214836589143940493046,
        8637771092812212464887027788957801177574860926032421582)


# it is useful if E's order is p.
def lift_up(point, a4, a6, p):
    # y^2=x^3+a4*x+a6 in Z/pZ
    E = EllipticCurve(Zmod(p ^ 2), [a4, a6])
    # \alpha(s,t)
    s, t = point.xy()
    # Z/p^2Z
    p2Z = Zmod(p ^ 2)
    # randomly choose a X1 mod p = s and it can be s it self
    X1 = p2Z(s)
    ###########generate Y1##########
    y = p2Z(t)
    omega = (((ZZ(X1) ^ 3+a4*ZZ(X1)+a6-ZZ(y) ^ 2)/p) % p)/(2*t)
    # if aseert work then pls change X1.
    assert omega != 0
    Y1 = y+p*p2Z(omega)
    ################################
    # X1, Y1 is a ZZ point in E(Z/p^2Z)
    A = E(X1, Y1)
    Xp_1, Yp_1 = ((p-1)*A).xy()
    # check if lambda_E is not zero
    assert X1 != Xp_1
    lambda_E = ((ZZ(Xp_1)-ZZ(X1))/p % p)/((ZZ(Yp_1)-ZZ(Y1)) % p)
    pZ = Zmod(p)
    # lambda_E is in Z/pZ
    lambda_E = pZ(lambda_E)
    return lambda_E


if __name__ == "__main__":
    lambda_P = lift_up(tP, A, B, p)
    lambda_Q = lift_up(tQ, A, B, p)
    c = lambda_Q/lambda_P
    if tQ == ZZ(c)*tP:
        print("ok!")
    print(f"secret_key: {c}")