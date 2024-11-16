
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from PIL import Image

import yaml
import pandas as pd

def loaddata(show_im=True):
    psf = Image.open(psfname)
    psf = np.array(psf, dtype='float32')
    data = Image.open(imgname)
    data = np.array(data, dtype='float32')

    """In the picamera, there is a non-trivial background 
    (even in the dark) that must be subtracted"""
    bg = np.mean(psf[5:15,5:15]) 
    psf -= bg
    data -= bg
    # capture psf
    #np.savetxt('real_psf_h_mem.txt', psf , fmt='%1.6f', delimiter= ' ')
    
    """Resize to a more manageable size to do reconstruction on. 
    Because resizing is downsampling, it is subject to aliasing 
    (artifacts produced by the periodic nature of sampling). Demosaicing is an attempt
    to account for/reduce the aliasing caused. In this application, we do the simplest
    possible demosaicing algorithm: smoothing/blurring the image with a box filter"""
    
    def resize(img, factor):
        num = int(-np.log2(factor))
        for i in range(num):
            img = 0.25*(img[::2,::2,...]+img[1::2,::2,...]+img[::2,1::2,...]+img[1::2,1::2,...])
        return img    
    
    psf = resize(psf, f)
    data = resize(data, f)

    """ To match the FFT that we will be using we crudely truncate the columns to match rows"""
    psf = psf[1:256, 1:256]
    data = data[1:256, 1:256]
    
    """ nmormalizing copy from shreyas"""
    psf /= np.linalg.norm(psf.ravel())
    data /= np.linalg.norm(data.ravel())
    
    if show_im:
        fig1 = plt.figure()
        plt.imshow(psf, cmap='gray')
        plt.title('PSF')
        plt.show()
        fig2 = plt.figure()
        plt.imshow(data, cmap='gray')
        plt.title('Raw data')
        plt.show()
    return psf, data

def initMatrices(h):
    pixel_start = (np.max(h) + np.min(h))/2
    x = np.ones(h.shape)*pixel_start

    init_shape = h.shape
    padded_shape = [nextPow2(2*n - 1) for n in init_shape]
    starti = (padded_shape[0]- init_shape[0])//2
    endi = starti + init_shape[0]
    startj = (padded_shape[1]//2) - (init_shape[1]//2)
    endj = startj + init_shape[1]
    hpad = np.zeros(padded_shape)
    hpad[starti:endi, startj:endj] = h

    H = fft.fft2(hpad, norm="ortho")
    Hadj = np.conj(H)

    def crop(X):
        cropped = X[starti:endi, startj:endj]
        return cropped

    def pad(v):
        vpad = np.zeros(padded_shape).astype(np.complex64)
        vpad[starti:endi, startj:endj] = v
        return vpad

    utils = [crop, pad]
    v = np.real(pad(x))
    
    return H, Hadj, v, utils

def nextPow2(n):
    return int(2**np.ceil(np.log2(n)))

def grad(Hadj, H, vk, b, crop, pad):
    Av = calcA(H, vk, crop)

    # Debug for determining bits carried
    # but clamping does not account for internal proc!!
    re_value, imag_value = split_complex(Av)
    re_clamp  =  clamp(re_value,1e-07,1,255)
    imag_clamp =  clamp(imag_value,1e-07,1,255)
    Av = re_clamp + imag_clamp

    diff = Av - b
    return np.real(calcAHerm(Hadj, diff, pad,crop))

def calcA(H, vk, crop):
    Vk = fft.fft2(vk, norm="ortho")

    #return crop(fft.ifftshift(fft.ifft2(Vk, norm="ortho")))

    # Debug for determining bits carried
    Had_prod = H*Vk
    ifft2_out = fft.ifftshift(fft.ifft2(Had_prod, norm="ortho"))
    return crop(ifft2_out)


def calcAHerm(Hadj, diff, pad,crop):
    xpad = pad(diff)
    X = fft.fft2(xpad, norm="ortho")

    #return fft.ifftshift(fft.ifft2(X, norm="ortho"))

    # Debug for determining bits carried
    Hadj_prod = Hadj*X
    #df = pd.DataFrame(Hadj_prod)
    #df.to_csv(r"C:\design\mig_zcu_104_2\fista_debug_data\input_to_ifft2.csv",index=False)
    #df  = pd.DataFrame(Hadj)
    #df.to_csv(r"C:\design\mig_zcu_104_2\fista_debug_data\hadj_data.csv",index=False)

    CHerm = fft.ifftshift(fft.ifft2(Hadj_prod, norm="ortho"))
    return CHerm


def clamp(n,minn,maxn,size):
    rows = size
    cols = size
    for i in range(rows):
        for j in range(cols):
            if n[i][j] < minn:
               n[i][j] = minn
            else:
                if n[i][j] > maxn:
                    n[i][j] = maxn


    return n

def split_complex(split_input):
    re_value = split_input.real
    imag_value = split_input.imag

    return re_value, imag_value

def grad_descent(h, b):
    H, Hadj, v, utils = initMatrices(h)
    crop = utils[0]
    pad = utils[1]
    # capture H
    np.savetxt('psf_big_H_mem.txt', H, fmt='%1.6f', delimiter=' ')

    alpha = np.real(2/(np.max(Hadj * H)))
    iterations = 0
     
    def non_neg(xi):
        xi = np.maximum(xi,0)
        return xi

    #proj = lambda x: x #Do no projection
    proj = non_neg #Enforce nonnegativity at every gradient step. Comment out as needed.


    parent_var = [H, Hadj, b, crop, pad, alpha, proj]
    
    vk = v
    
    
    
    #### uncomment for Nesterov momentum update ####   
    #p = 0
    #mu = 0.9
    ################################################
    
    
    
    #### uncomment for FISTA update ################
    tk = 1
    xk = v
    ################################################
        
    for iterations in range(iters):   
        
        # uncomment for regular GD update
        #vk = gd_update(vk, parent_var)
        
        # uncomment for Nesterov momentum update 
        #vk, p = nesterov_update(vk, p, mu, parent_var)
        
        # uncomment for FISTA update
        vk, tk, xk = fista_update(vk, tk, xk, parent_var)

        if iterations % disp_pic == 0:
            print(iterations)
            image = proj(crop(vk))
            f = plt.figure(1)
            plt.imshow(image, cmap='gray')
            plt.title('Reconstruction after iteration {}'.format(iterations))
            plt.show()
    
    
    return proj(crop(vk)) 
    
def gd_update(vk, parent_var):
    H, Hadj, b, crop, pad, alpha, proj = parent_var
    
    gradient = grad(Hadj, H, vk, b, crop, pad)
    vk -= alpha*gradient
    vk = proj(vk)
    
    return xk    

def nesterov_update(vk, p, mu, parent_var):
    H, Hadj, b, crop, pad, alpha, proj = parent_var
    
    p_prev = p
    gradient = grad(Hadj, H, vk, b, crop, pad)
    p = mu*p - alpha*gradient
    vk += -mu*p_prev + (1+mu)*p
    vk = proj(vk)
    
    return vk, p

def fista_update(vk, tk, xk, parent_var):
    H, Hadj, b, crop, pad, alpha, proj = parent_var
    
    x_k1 = xk
    gradient = grad(Hadj, H, vk, b, crop, pad)
    vk -= alpha*gradient
    #xk = proj(vk)
    #t_k1 = (1+np.sqrt(1+4*tk**2))/2
    #vk = xk+(tk-1)/t_k1*(xk - x_k1)
    #tk = t_k1
    
    return vk, tk, xk


if __name__ == "__main__":
    ### Reading in params from config file (don't mess with parameter names!)
    #params = yaml.load(open("gd_config.yml"))
    #for k,v in params.items():
    #    exec(k + "=v")
    psfname = "./test_images/PSF_RULER_PINHOLE.tif"
    imgname = "./test_images/ThorlabsDog_Image__2021-12-23__11-07-19.tif"
    #psfname = "./test_images/psf_sample.tif"  # Path to PSF image
    #imgname = "./test_images/rawdata_hand_sample.tif"  # Path to raw data image
    #psfname = "./test_images/double_side_xy_psf.tif"  # Path to PSF image
    #imgname = "./test_images/double_side_xy_image.tif"  # Path to raw data image
    #psfname = "./test_images/single_side_xy_psf.tif"  # Path to PSF image
    #imgname = "./test_images/single_side_xy_image.tif"  # Path to raw data image
    f = 0.125  # Downsampling factor (must be decimal, must be 1/2^k where k is positive integer)
    iters = 120  # Number of iterations
    #iters = 20
    disp_pic = 20  # Number of iterations after which we display intermediate reconstruction

    psf, data = loaddata()
    final_im = grad_descent(psf, data)
    print(iters)
    plt.imshow(final_im, cmap='gray')
    plt.title('Final reconstruction after {} iterations'.format(iters))
    plt.show()
    saveim = input('Save final image? (y/n) ')
    if saveim == 'y':
        filename = input('Name of file: ')
        plt.imshow(final_im, cmap='gray')
        plt.axis('off')
        plt.savefig(filename+'.png', bbox_inches='tight')

