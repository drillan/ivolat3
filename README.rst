Metadata-Version: 1.1
Name: ivolat3
Version: 0.0.0
Author: driller
Home-page: https://github.com/drillan/ivolat3
Summary: European Options Pricing Library
License: MIT
Description: Installation
        ------------
        ::
        
            pip install ivolat3
        
        
        Examples
        --------
        ::
        
            import ivolat3
            
            # stock price
            s = 10000
            # option strike price
            k = 10500
            # risk-free interest rate
            r = 0.001
            # drift rate (dividend)
            q = 0.1
            # time remaining until expiration (in years)
            t = 20.0 / 365.0
            # annual volatility of stock price
            sigma = 0.2
            # call option price
            p = 40
        
            
            # call opition price
            ivolat3.prem_call(s, k, r, q, t, sigma)
            
            # call opition implied volatility
            ivolat3.ivolat_call(s, k, r, q, t, p)
            
            # Greeks
            # call delta
            ivolat3.delta_call(s, k, r, q, t, sigma)
            
            # Second-order Greeks
            # gamma
            ivolat3.gamma(s, k, r, q, t, sigma)
            
            # vanna
            ivolat3.vanna(s, k, r, q, t, sigma)
            
            # vomma (volga)
            ivolat3.vomma(s, k, r, q, t, sigma)
Platform: UNKNOWN
Classifier: Operating System :: Microsoft :: Windows
Classifier: Operating System :: POSIX :: Linux
Classifier: Programming Language :: Python :: 3.5