try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

VERSION = "0.1.2"
LONG_DESCRIPTION = """
calculate items by Black-Scholes model.
* call price
* put price
* implied volatility
* delta
* vega
* theta
* gamma
* vanna
* charm
* speed
* zomma
* color
* DvegaDtime
* vomma
* ultima
* dualdelta
* dualgamma
"""

if __name__ == '__main__':
    bs = Extension('bs', sources=['ivolat3/bs.c'])
    setup(name="ivolat3", version=VERSION,
          url='https://github.com/drillan/ivolat3',
          description="European Options Pricing Library",
          long_description=LONG_DESCRIPTION,
          classifiers=[
              'Programming Language :: Python :: 3',
          ],
          packages=['ivolat3'], ext_modules=[bs])
