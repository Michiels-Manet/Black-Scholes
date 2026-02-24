using System; // console.writeline ipv system.console.writeline

public static class BlackScholes // static = geen object aanmaken
{
    // Europese call via Black–Scholes
    public static double CallPrice(double S, double K, double r, double sigma, double T) // S = huidige prijs van het onderliggende actief, K = uitoefenprijs, r = risicovrije rente, sigma = volatiliteit, T = tijd tot maturity
    {
        ValidateInputs(S, K, sigma, T); // controleert of S,K,sigma,T positief zijn, anders error

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
        double d2 = d1 - sigma * sqrtT;

        return S * NormCdf(d1) - K * Math.Exp(-r * T) * NormCdf(d2); // NormCdf = cummulatieve normale verdeling
    }

    // Europese put via Black–Scholes vergelijking
    public static double PutPrice(double S, double K, double r, double sigma, double T)
    {
        ValidateInputs(S, K, sigma, T);

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
        double d2 = d1 - sigma * sqrtT;

        return K * Math.Exp(-r * T) * NormCdf(-d2) - S * NormCdf(-d1);
    }

    private static void ValidateInputs(double S, double K, double sigma, double T)
    {
        if (S <= 0) throw new ArgumentOutOfRangeException(nameof(S), "S moet > 0 zijn."); // log(S/K) bestaat niet
        if (K <= 0) throw new ArgumentOutOfRangeException(nameof(K), "K moet > 0 zijn."); // log(S/K) bestaat niet en delen door 0 gaat niet
        if (sigma <= 0) throw new ArgumentOutOfRangeException(nameof(sigma), "sigma moet > 0 zijn."); // delen door 0 gaat niet
        if (T <= 0) throw new ArgumentOutOfRangeException(nameof(T), "T moet > 0 zijn (in jaren)."); // sqrt(T) bestaat niet
    } 

    // Standaard normale CDF Φ(x) zonder Erf:
    // veelgebruikte numerieke benadering (Abramowitz/Stegun-achtige polynomiale benadering)
    private static double NormCdf(double x)
    {
        double absX = Math.Abs(x);
        double t = 1.0 / (1.0 + 0.2316419 * absX);

        double d = 0.3989422804014327 * Math.Exp(-0.5 * absX * absX);

        double prob = d * t * (
            0.319381530 +
            t * (-0.356563782 +
            t * (1.781477937 +
            t * (-1.821255978 +
            t * 1.330274429)))
        );

        return (x >= 0) ? (1.0 - prob) : prob; // symmetrie van de normale verdeling
    }
}

class Program
{
    static void Main(string[] args) // hier begint je programma te runnen
    {

        // Testwaarden - voorbeeldparameters
        double S = 100.0;
        double K = 100.0;
        double r = 0.05;
        double sigma = 0.20;
        double T = 1.0;

        double call = BlackScholes.CallPrice(S, K, r, sigma, T); // bereken call
        double put = BlackScholes.PutPrice(S, K, r, sigma, T); // bereken put

        // print de resultaten
        Console.WriteLine($"Call prijs: {call:F6}"); // F6 = 6 decimalen
        Console.WriteLine($"Put  prijs: {put:F6}");
    }
}
