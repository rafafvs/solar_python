from .boundTransform import BoundTransform

class SolarTransform(BoundTransform):
    """Solar Model transformation functions connecting Clear Sky to Yt."""
    
    def __init__(self, alpha=0.0, beta=1.0, link="invgumbel"):
        super().__init__(alpha, beta, link)

    def X(self, Rt, Ct):
        """Map the solar radiation Rt in the risk driver Xt."""
        return 1.0 - (Rt / Ct)

    def iX(self, Xt, Ct):
        """Map the risk driver Xt in solar radiation Rt."""
        return Ct * (1.0 - Xt)

    def eta(self, Rt, Ct):
        """Map the solar radiation Rt in the normalized variable Xt_prime."""
        return self.X_prime(self.X(Rt, Ct))

    def ieta(self, Xt_prime, Ct):
        """Map the normalized variable Xt_prime to the solar radiation Rt."""
        return self.iX(self.iX_prime(Xt_prime), Ct)

    def RY(self, Rt, Ct):
        """Convert solar radiation Rt into the transformed variable Yt."""
        Xt = self.X(Rt, Ct)
        Xt_prime = self.X_prime(Xt)
        return self.Y(Xt_prime)

    def iRY(self, Yt, Ct):
        """Convert the transformed variable Yt into solar radiation Rt."""
        Xt_prime = self.iY(Yt)
        Xt = self.iX_prime(Xt_prime)
        return self.iX(Xt, Ct)

    def __str__(self):
        """Print method for the class `SolarTransform`."""
        msg = f"------------------------ \033[1;35m boundTransform \033[0m ------------------------", "\n"
        msg += f" - Link: \033[1;32m {self.link} \033[0m \n"
        msg += f" - iY(y): \033[1;32m {self.alpha:.3f} + \033[1;32m {self.beta:.4f} exp(-exp(Y)) \033[0m \n"
        msg += f" - Y(x): log(log(\033[1;32m {self.beta:.4f} \033[0m) - log(x- \033[1;32m {self.alpha:.3f} \033[0m))"

        return msg

    def __repr__(self):
        return self.__str__()