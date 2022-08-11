p = 730750818665451459112596905638433048232067471723
A = 425706413842211054102700238164133538302169176474
B = 203362936548826936673264444982866339953265530166

# y^2 = x^3 + A * x + B
tE = EllipticCurve(GF(p), [A, B])

tP = tE(125270202464411072778547771568975423382990845665,
        440970603958123875213441435758390311809187352362)
tQ = tE(296939092187233862778999244256460019221379646447,
        650021996391906816753000782344884033217784284632)


# it is useful if E's order is p.
def lift_up(point, tE, p):
    # y^2=x^3+a4*x+a6 in Z/pZ
    E = tE.change_ring(Zp(p, 2))
    a4 = tE.a4()
    a6 = tE.a6()
    # \alpha(s,t)
    s, t = point.xy()

    # Zp
    zp = Zp(p, 3)
    # randomly choose a X1 mod p = s and it can be s it self
    X1 = zp(s)
    ###########generate Y1##########
    y = zp(t)
    omega = (((ZZ(X1) ^ 3+ZZ(a4)*ZZ(X1)+ZZ(a6)-ZZ(y) ^ 2)/p) % p)/(2*ZZ(t))
    # if aseert work then pls change X1.
    assert omega != 0
    Y1 = y+p*zp(omega)
    ################################
    # X1, Y1 is a  point in E(Z/p^2Z)
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
    lambda_P = lift_up(tP, tE, p)
    lambda_Q = lift_up(tQ, tE, p)
    c = lambda_Q/lambda_P
    if tQ == (c)*tP:
        print("ok!")
    print(f"secret_key: {c}")
