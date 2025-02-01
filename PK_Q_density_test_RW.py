# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:07:59 2024

@author: leizhou
"""

def nonneg_pseudoinverse(A):

     m, n = A.shape
     I = np.eye(n)  # 单位矩阵
     X = np.zeros((n, n))  # 初始化结果矩阵
     
     for j in range(n):  
         X[:, j], _ = nnls(A, I[:, j])  # 求解 Ax = b, 使得 x >= 0
     
     return X
 
def Ker(u):
    # Gaussian kernel function (standard normal distribution)
    if isinstance(u, (int, float)):  # If u is a scalar
        return norm.pdf(u)  # Calculate PDF for scalar u
    else:
        return np.array([norm.pdf(val) for val in u])  # Handle list or array of u

def K_h(u, h):
    # Scaled kernel function: K(u/h) / h
    if isinstance(u, (int, float)):  # If u is scalar
        return Ker(u / h) / h
    else:
        return np.array([Ker(val / h) / h for val in u])  # Handle list or array of u

def K_prime(u):
    # Derivative of the Gaussian kernel
    if isinstance(u, (int, float)):  # If u is scalar
        return -u * norm.pdf(u)  # Calculate derivative for scalar u
    else:
        return np.array([-val * norm.pdf(val) for val in u])  # Handle list or array of u

def Kjp(u, j, N):
    # Kernel function with polynomial terms and determinant adjustment
    kk = len(u)
    ret = np.zeros(kk)
    
    for i in range(kk):
        M = N.copy()
        M[:, j] = np.power(u[i], np.arange(N.shape[0]))
        
        if abs(u[i]) > 1:
            ret[i] = 0
        else:
            det_M = np.linalg.det(M)  # Determinant of matrix M
            det_N = np.linalg.det(N)  # Determinant of matrix N
            ret[i] = Ker(u[i]) * det_M / det_N if det_N != 0 else 0
            
    return ret

def func_L(p, h, alpha, n):  # Add `n` as a parameter for sample size or related value
    
    # Initialize matrices and vectors
    N = np.zeros((p+1, p+1))
    Q = np.zeros((p+1, p+1))
    C = np.zeros(p+1)
    L = np.zeros(p+1)

    # Helper function for integration
    def integrate_func(func, lower, upper):
        result, _ = quad(func, lower, upper)
        return result

    # Fill matrix Q and N with appropriate integrals
    for i in range(p + 1):  # Change to range(p + 1) to include p
        for j in range(p + 1):  # Change to range(p + 1) to include p
            if i > 1 and j > 1:
                Q[i, j] = integrate_func(lambda xx: xx**(i + j) * (K_prime(xx))**2, -np.inf, np.inf) - 0.5 * (i * (i - 1) + j * (j - 1)) * integrate_func(lambda x: x**(i + j - 2) * (Ker(x))**2, -np.inf, np.inf)
            elif i <= 1 and j <= 1:
                Q[i, j] = integrate_func(lambda xx: xx**(i + j) * (K_prime(xx))**2, -np.inf, np.inf)
            elif i <= 1 and j > 1:
                Q[i, j] = integrate_func(lambda xx: xx**(i + j) * (K_prime(xx))**2, -np.inf, np.inf) - 0.5 * j * (j - 1) * integrate_func(lambda x: x**(i + j - 2) * (Ker(x))**2, -np.inf, np.inf)
            elif i > 1 and j <= 1:
                Q[i, j] = integrate_func(lambda xx: xx**(i + j) * (K_prime(xx))**2, -np.inf, np.inf) - 0.5 * i * (i - 1) * integrate_func(lambda x: x**(i + j - 2) * (Ker(x))**2, -np.inf, np.inf)
            N[i, j] = integrate_func(lambda xx: xx**(i + j) * Ker(xx), -1, 1)


    # Invert matrix N
    N_i = np.linalg.inv(N)

    # Compute C and L
    for j in range(p + 1):
        C[j] = (N_i @ Q @ N_i)[j, j] / integrate_func(lambda xx: (Kjp([xx], j, N))**2, -np.inf, np.inf)
        L[j] = factorial(j) * (n * h**(2 * j + 1))**(-0.5) * ((-2 * log(h))**0.5 + (-2 * log(h))**(-0.5) * (-log(-0.5 * log(1 - alpha)) + log(sqrt(C[j]) / (2 * pi))))

    print(C)
    return L

def bands(p, h, hat_theta, x, y, vx):
    # p=2
    # h = 0.08
    # x =np.array(sub['moneyness']) 
    # y = np.array(sub['mark_iv']) 
    # vx = M_std
    # hat_theta = sigmas
    
    L = func_L(p, h, 0.05, n = len(x))
    sigma2=1
       
    up_band = np.zeros((len(vx), p + 1))
    lo_band = np.zeros((len(vx), p + 1))
    V_all = np.zeros((len(vx), p + 1))
    
    for i in range(len(vx)):
        xx = vx[i]
        B_V = np.zeros((p + 1, p + 1))
        K_V = np.zeros((p + 1, p + 1))
    
        for j in range(len(x)):
            # 计算 K.h(x[j] - xx, h)
         
            K_h_val = K_h(x[j] - xx, h)
    
            # 计算对角矩阵 diag(h^(0:p))
            diag_matrix = np.diag(h ** np.arange(0, p + 1))
            inv_diag_matrix = np.linalg.inv(diag_matrix)
    
            # 计算 (x[j] - xx) ** (0:p) 并确保它是列向量
            x_diff_power = np.power(x[j] - xx, np.arange(0, p + 1)).reshape(-1, 1)  # 列向量
    
            # 更新 B.V
            B_V += K_h_val * (1 / sigma2) * (inv_diag_matrix @ (x_diff_power @ x_diff_power.T) @ inv_diag_matrix)
    
            # 计算 K.V
            term1 = K_h_val ** 2
            residual = y[j] - np.sum(hat_theta[i, :] * (x[j] - xx) ** np.arange(p + 1))
            term2 = (residual / sigma2) ** 2
    
            # 计算矩阵
            mat_product = inv_diag_matrix @ (x_diff_power @ x_diff_power.T) @ inv_diag_matrix
    
            # 更新 K.V
            K_V += h * term1 * term2 * mat_product
    
        # 归一化
        B_V /= len(x)
        K_V /= len(x)
    
        # 计算 V
        B_inv = np.linalg.inv(B_V)
        V = B_inv @ K_V @ B_inv
        V_all[i, :] = np.diag(V)
    
        # 计算置信区间
        for ip in range(p + 1):
            lo_band[i, ip] = -L[ip] * np.sqrt(V[ip, ip])
            up_band[i, ip] = L[ip] * np.sqrt(V[ip, ip])

    return {'L': L, 'V': V_all, 'lo_band': lo_band, 'up_band': up_band}

def create_single_folder(folder_path):
    try:
        os.mkdir(folder_path)
        print(f"file '{folder_path}' create")
    except FileExistsError:
        print(f" '{folder_path}' exists")
    except Exception as e:
        print(f" '{folder_path}' wrong: {e}")



def gaussian_kernel(M, m, h_m, T, t, h_t):
    u_m = (M-m)/h_m
    u_t = (T-t)/h_t
    return stats.norm.cdf(u_m) * stats.norm.cdf(u_t)

def epanechnikov(M, m, h_m, T, t, h_t):
    u_m = (M-m)/h_m
    u_t = (T-t)/h_t
    return (3/4) * (1-u_m)**2 * (3/4) * (1-u_t)**2

def extend_polynomial(x, y):
    """
    Extend Smile, first and second derivative so that spd exists completely for large tau
    x = M_std
    y = first
    """
    polynomial_coeff=np.polyfit(x,y,2)
    xnew=np.linspace(0.8,1.2,100)
    ynew=np.poly1d(polynomial_coeff)
    
    # plt.plot(xnew,ynew(xnew),x,y,'o')
    # plt.title('interpolated smile')
    # plt.show()
    return xnew, ynew(xnew)

def smoothing_rookley(df, m, t, h_m, h_t, IV_nn, M_nn):
    
    # kernel=gaussian_kernel, 
    extend = True
    
    # Before
    # M = np.array(df.moneyness)
    # y = np.array(df.mark_iv)
    
    M = M_nn
    y = IV_nn

    # After Poly extension
    if extend:
        print('Extending Moneyness and IV in smoothing technique!')
        M, y = extend_polynomial(M, y)
     
    
    # temp ={'IV':y, 'M':M}
    # temp = pd.DataFrame(temp)
    # plt.figure(figsize=(10, 6))
    # plt.scatter(temp['M'],temp['IV'], color='r')
    # plt.show   
     
    T = [df.tau.values[0]] * len(M) #np.array(df.tau)
    n = len(M)

    X1 = np.ones(n)
    X2 = M - m
    X3 = (M-m)**2
    X4 = T-t
    X5 = (T-t)**2
    X6 = X2*X4
    X = np.array([X1, X2, X3, X4, X5, X6]).T

    # ker = kernel(M, m, h_m, T, t, h_t)
    ker = gaussian_kernel(M, m, h_m, T, t, h_t)
    W = np.diag(ker)

    XTW = np.dot(X.T, W)

    beta = np.linalg.pinv(np.dot(XTW, X)).dot(XTW).dot(y)

    return beta[0], beta[1], 2*beta[2]


def rookley_unique_tau(df, h_m):
    # gridsize is len of estimated smile
    # df = sub
    # h_m = h
    
    h_t=0.01
    gridsize=149
    kernel='epak'

    if kernel=='epak':
        kernel = epanechnikov
    elif kernel=='gauss':
        kernel = gaussian_kernel
    else:
        print('kernel not know, use epanechnikov')
        kernel = epanechnikov

    
    num = gridsize
    M_min, M_max = min(df.moneyness), max(df.moneyness)
    # M = np.linspace(M_min, M_max, gridsize)
    
   
    # M_n = np.unique(df['moneyness'])
    # IV_n = np.full(len(M_n),np.nan)

    # for i, m in enumerate(M_n):
    #     mask = df['moneyness'] == m
    #     IV_n[i] = np.mean(df['mark_iv'][mask])
        
    M_nn = np.linspace(M_min, M_max, 1001)
    IV_nn = np.full(len(M_nn) - 1, np.nan)  # Make IV_nn one element shorter

    for i in range(len(M_nn) - 1):  # Iterate only up to the second last element
        mask = (df['moneyness'] >= M_nn[i]) & (df['moneyness'] < M_nn[i + 1])
        IV_nn[i] = np.mean(df['mark_iv'][mask])

    M_nn = M_nn[:-1]    
    
    M_nn, IV_nn = M_nn[~np.isnan(IV_nn)],IV_nn[~np.isnan(IV_nn)]
    

    # temp ={'IV':IV_nn, 'M':M_nn}
    # temp = pd.DataFrame(temp)
    # plt.figure(figsize=(10, 6))
    # plt.scatter(temp['M'],temp['IV'], color='r')
    # plt.show
 
    
    # M_std_min, M_std_max = min(df.moneyness), max(df.moneyness)
    # M_std = np.linspace(M_std_min, M_std_max, num=num)
    
    # M_std_min, M_std_max = min(M_nn), max(M_nn)
    M_std_min, M_std_max = 0.8, 1.2
    M_std = np.linspace(M_std_min, M_std_max, num=num)

    # if all taus are the same
    tau_min, tau_max = min(df.tau[(df.tau > 0)]), max(df.tau) # empty sequence for tau precision = 3
    tau = np.linspace(tau_min, tau_max, gridsize)

    x = zip(M_std, tau)
    sig = np.zeros((num, 3)) # fill
    

    # TODO: speed up with tensor instead of loop
    # (m, t) here means (M_std and tau) in x
    for i, (m, t) in enumerate(x):
        sig[i] = smoothing_rookley(df, m, t, h_m, h_t, IV_nn, M_nn)
  
    
    smile = sig[:, 0]
    first = sig[:, 1] #/ np.std(df.moneyness)
    second = sig[:, 2] #/ np.std(df.moneyness)

    # S_min, S_max = min(df.index_price), max(df.index_price)
    # K_min, K_max = min(df.strike), max(df.strike)
    # S = np.linspace(S_min, S_max, gridsize)
    # K = np.linspace(K_min, K_max, gridsize)
    S = [df['index_price'].mean() for _ in range(len(M_std))]
    K = S / M_std

    return smile, first, second, S, K, M_std, tau   

def compute_spd(sigma, sigma1, sigma2, tau, s, m2, r_const):
    
    # tau = newtau
    # s = S
    # m2 = M_std

    # SPDBL
    # Scalars
    #tau = np.mean(tau)
    #s = np.mean(s)
    r = r_const * len(s)
    #r = np.mean(r)
    #r = 0

    # now start spdbl estimation
    st = np.sqrt(tau)
    ert = np.exp(r * tau)
    rt = r * tau
    # error should be here in the length of m
    d1 = (np.log(m2) + tau * (r + 0.5 * (sigma ** 2))) / (sigma * st)
    d2 = d1 - (sigma * st)
    f = stats.norm.cdf(d1, 0, 1) - (stats.norm.cdf(d2, 0, 1)/(ert * m2))

    # First derivative of d1
    d11 = (1/(m2*sigma*st)) - (1/(st*(sigma**2))) * ((np.log(m2) + tau * r) * sigma1) + 0.5 * st * sigma1
    
    # First derivative of d2 term
    d21 = d11 - (st * sigma1)
    
    # Second derivative of d1 term
    d12 = -(1/(st * (m2**2) * sigma)) - sigma1/(st * m2 * (sigma**2)) + sigma2 * (0.5 * st - (np.log(m2) + rt)) / (st * sigma**2) + sigma1 * (2 * sigma1 * (np.log(m2) + rt)) / (st * sigma**3) - 1/(st * m2 * sigma**2)

    # Second derivative of d2 term
    d22 = d12 - (st * sigma2)

    f1 = (stats.norm.pdf(d1, 0, 1) * d11) + (1/ert) * ((-stats.norm.pdf(d2, 0, 1) * d21)/m2 + stats.norm.cdf(d2, 0, 1) / m2**2)
    
    # f2 = dnorm(d1, mean = 0, sd = 1) * d12 - d1 * dnorm(d1, mean = 0, sd = 1) * (d11^2) - (1/(ert * m) * dnorm(d2, mean = 0, sd = 1) * d22) + ((dnorm(d2, mean = 0, sd = 1) * d21)/(ert * m^2)) + (1/(ert * m) * d2 * dnorm(d2, mean = 0, sd = 1) * (d21^2)) - (2 * pnorm(d2, mean = 0, sd = 1)/(ert * (m^3))) + (1/(ert * (m^2)) * dnorm(d2, mean = 0, sd = 1) * d21)
    f2 = stats.norm.pdf(d1, 0, 1) * d12 - d1 * stats.norm.pdf(d1, 0, 1) * d11**2     - (1/(ert * m2) * stats.norm.pdf(d2, 0, 1) * d22) + ((stats.norm.pdf(d2, 0, 1)*d21)/(ert * m2**2))      + (1/(ert * m2) * d2 * stats.norm.pdf(d2, 0, 1)) * d21**2 -(2 * stats.norm.cdf(d2, 0, 1)/(ert * m2**3))     + (1/(ert * m2**2)) * stats.norm.pdf(d2, 0, 1) *d21
   
    # recover strike price
    x = s/m2
    c1 = -(m2**2) * f1
    c2 = s * (1/x**2) * ((m2**2) * f2 + 2 * m2 * f1)

    # Calculate the quantities of interest
    cdf = ert * c1 + 1
    fstar = ert * c2
    delta = f + s + f1/x
    gamma = 2 * f1 / x + s * f2 / (x**2)
    #print('\ndelta: ', delta, '\ngamma:', gamma)

    #plt.plot(fstar)
    #plt.show()

    spd = pd.DataFrame({'x': np.mean(s)/m2, # strike
                        'y': fstar, # Q density
                        'm': m2}) # Moneyness


    return spd  


def spdbl(sub, date, tau, r_const):
    
    # spd, sub = spdbl(sub_pre, date, tau, int_rate)
    
 
    # r_const = int_rate
    """
    computes spd according to Breeden and Litzenberger 1978
    """
    # plot_ident = '2020-' + str(mindate.month) + '-' + str(mindate.day) + '-' + str(tau) + '.png'
    
    tau_file = tau
    
    # Subset
    # Only calls, tau in [0, 0.25] and fix one day (bc looking at intra day here)
    
  
        # 
        # @Todo: tau should always be > 0, else check!
    
    if sub.shape[0] == 0:
        raise(ValueError('Sub is empty'))

    del sub['date']
    #sub['tau'] = round(sub['tau'], 2)
    sub['moneyness'] = round(sub['moneyness'], 3)
    sub['index_price'] = round(sub['index_price'], 2)

    sub = sub.drop_duplicates()
    #print('Only unique transactions!')
    
    print(sub.describe())

    if sub.shape[1] > 45000:
        raise(ValueError('Sub too large'))

    # Isolate vars
    sub['mark_iv'] = sub['mark_iv']/100
    # sub['mark_iv'][(sub['mark_iv'] < 0.01)] = 0
    
    M_range =  np.arange(sub['moneyness'].min(),sub['moneyness'].max() , 0.05)
   
    nrows = 2
    ncols = ceil(len(M_range) / nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 6 * nrows))
    
    # 将axes转换为1D数组，方便迭代
    axes = axes.flatten()
    
    # 记录有效子图数量
    num_plots = 0
    
    sub_keep = pd.DataFrame()
    Sum_IV = pd.DataFrame()
    for i, M_t in enumerate(M_range):
        Scale = sub[(sub['moneyness'] >= M_t) & (sub['moneyness'] < M_t + 0.05)]
        if Scale.empty:
            continue  # 如果数据为空，跳过这个循环
            
        mean_iv = Scale['mark_iv'].mean()
        mean_M = Scale['moneyness'].mean()
       
        lower_quantile = Scale['mark_iv'].quantile(0.03)
        median_quantile = Scale['mark_iv'].quantile(0.5)
        upper_quantile = Scale['mark_iv'].quantile(0.97)
    
        sns.boxplot(y='mark_iv', data=Scale, ax=axes[num_plots])
        axes[num_plots].set_ylabel('Values')
        num_observations = len(Scale)
        
        out = pd.DataFrame({
                'M_start':[Scale['moneyness'].min()],
                'M_end':[Scale['moneyness'].max()],
               'mean_M': [mean_M], 
               'obs': [num_observations], 
               'mean_iv': [mean_iv], 
               'sd_IV': [Scale['mark_iv'].std()], 
               'min_IV': [Scale['mark_iv'].min()], 
               'max_IV': [Scale['mark_iv'].max()]
               })
       
        Sum_IV =  pd.concat([Sum_IV,out])
        
        
        axes[num_plots].set_title(f'Moneyness {round(M_t, 3)} to {round(M_t + 0.05, 3)}\nN={num_observations}')
    
        # 添加分位线
        axes[num_plots].axhline(y=lower_quantile, color='r', linestyle='--', label='3rd Percentile')
        axes[num_plots].axhline(y=median_quantile, color='g', linestyle='--', label='Median (50th Percentile)')
        axes[num_plots].axhline(y=upper_quantile, color='b', linestyle='--', label='97th Percentile')
        
        Scale_keep = Scale[(Scale['mark_iv']<upper_quantile) & (Scale['mark_iv']>lower_quantile)]
        
        if Scale_keep.empty == 0:
            sub_keep = pd.concat([sub_keep,Scale_keep])
        
        num_plots += 1
    
    # 删除未使用的子图
    for j in range(num_plots, len(axes)):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(path+'/Box_plot/Crypto_option/Box_'+np.datetime_as_string(date)[:10]+'_tau_'+str(int(tau_file))+'.png', transparent=True)
    
    # 显示图形
    plt.show()
    plt.close()
    
  
    
    temp ={'IV':sub['mark_iv'], 'M':sub['moneyness']}
    temp = pd.DataFrame(temp)
    plt.figure(figsize=(10, 6))
    plt.scatter(temp['M'],temp['IV'], color='r')
    plt.show
    plt.close()
    
    output_filename =  path+'/Box_plot/Crypto_summary/IV_'+np.datetime_as_string(date)[:10]+'_tau_'+str(int(tau_file))+'.csv'
    Sum_IV.to_csv(output_filename, index=False)
   
    
    # omit strange IV
    # sub = sub[(sub['mark_iv'] > lower_quantile) & (sub['mark_iv'] < upper_quantile)]
    sub = sub_keep
   
    
    if (sub.shape[0] == 0) | (np.unique(sub['moneyness']).shape[0] < 5):
        spd = 0
        return spd,sub
    else:
        
        h = sub.shape[0] ** (-1 / 9)
    
        # vola = sub['mark_iv'].astype(float)/100
        tau = sub['tau'].astype(float)
        tau_uni = np.array(tau)[0]
        #m = float(sub['index_price']/sub['strike'] # Spot price corrected for div divided by k ((s-d)/k); 
        #moneyness of an option; div always 0 here
        r = sub['interest_rate'].astype(float)
        s = sub['index_price'].astype(float)
        k = sub['strike'].astype(float)
        m = s / k
        div = 0
    
        
        # Forward price
        F = s * np.exp((r-div) * tau)
        K = m * F # k capital
        newm =(s * np.exp(-div*tau))/K
    
        sigma = []
        sigma1 = []
        sigma2 = []
    
        
        print('Choosing Bandwidth h: ', h)
        
        sigma, sigma1, sigma2, S, K, M_std, newtau = rookley_unique_tau(sub, h) # rookley(sub, h) 
    
        # save data
        temp ={'sigma':sigma, 'sigma1':sigma1, 'sigma2':sigma2,'S': S, 'K':K, 'M_std':M_std, 'newtau':newtau}
        temp = pd.DataFrame(temp)
        plt.figure(figsize=(10, 6))
        plt.scatter(temp['M_std'],temp['sigma'], color='r')
        plt.show
        
        plt.savefig(path+'Rookley_IV/'+np.datetime_as_string(date)[:10]+'_tau_'+str(int(tau_file))+'.png', transparent=True) 
        plt.close('all')  
        # sub2 = sub[(sub['mark_iv']>0.5) & (sub['mark_iv']<0.8)]
        # plt.scatter(sub2['moneyness'],sub2['mark_iv'], color='r')
        
        
        # tau is too long here!!
        spd = compute_spd(sigma, sigma1, sigma2, newtau,  S, M_std, r_const)
        
        sigmas = [sigma, sigma1, sigma2]
        sigmas = np.array(sigmas)
        sigmas = sigmas.T
        
        bands_init = bands(p=2, h=0.08, hat_theta=sigmas, x=np.array(sub['moneyness']) , y=np.array(sub['mark_iv']) , vx=M_std)
       
        # d1=(log(mon.grid.scaled)+tau*(int_rate+0.5*(sigma)^2)) / (sqrt(sigma)*sqrt(tau));
        d1 = (np.log(M_std) + tau_uni * (int_rate + 0.5 * (sigma ** 2))) / (sigma * np.sqrt(tau_uni))
        
        d2=d1-np.sqrt(sigma)*sqrt(tau_uni);
  
        dgds = np.exp(-int_rate * tau_uni) * (M_std ** 2 * norm.pdf(d1) / spd['x'] ** 2 - np.exp(-int_rate * tau_uni) * norm.pdf(d2) / M_std)
       
        band_limit = bands_init["L"][2] * np.sqrt(bands_init["V"][:,2]) * dgds
        band_lo = spd['y'] - abs(spd['y']*band_limit)
        band_up = spd['y'] + abs(spd['y']*band_limit)
        spd_band = pd.DataFrame({
                'm': spd['m'],  
                'Strike': spd['x'],  
                'SPD': spd['y'],    
                'SPD_lo': band_lo,  
                'SPD_up': band_up   
            })
        
        # plt.figure(figsize=(10, 6))
        # plt.scatter(spd['m'],spd['y'], color='r')
        # plt.show
        
        
        # kde = gaussian_kde(spd['m'], weights=spd['y'])
        # M_kde = np.linspace(max(spd['m']), 1.2, 100)
        # Q_kde = kde(M_kde)
        
        # degree = 4
        # coefficients = np.polyfit(spd['m'], spd['y'], degree)
        # polynomial = np.poly1d(coefficients)

        # M_fit = np.linspace(max(spd['m']), 1.2, 1001)
        # Q_fit = polynomial(M_fit)
        
        # # 生成拟合曲线的值
  
        # # 补充1.05到1.2区间的数据
        # M_combined = np.concatenate((spd['m'], M_fit))
        # Q_combined = np.concatenate((spd['y'], Q_fit))

        # plt.figure(figsize=(10, 6))
        # plt.scatter(M_combined,Q_combined, color='r')
        # plt.show
    
    return spd, sub, spd_band


# main body

import os
import pandas as pd
# import datetime
# import statsmodels.api as sm
# import matplotlib.pyplot as plt
import numpy as np 
# import time
# import sys
# import pdb
# import pickle
# import gc
# import csv
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil


from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from math import factorial, pi, log, sqrt
from scipy.stats import norm
from scipy.optimize import nnls

# path = '/Users/tracy/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/PAPER/Pricing/code/test/'
path = "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/main_code/"
os.chdir(path)

create_single_folder(path+'Rookley_IV')
# file = 'final_btc_option_2021-01-01-2022-03-31.csv'
# file = 'final_btc_option_2021-04-01-2021-06-30.csv'
# file = 'final_btc_option_2021-07-01-2021-09-29.csv'
# file = 'final_btc_option_2021-10-01-2021-12-30.csv'
# file = 'final_btc_option_2022-07-01-2022-09-29.csv'
file = 'final_btc_option_2022-10-01-2022-12-30.csv'
filename = os.path.join( path,'BRC_final/', file)
data_option = pd.read_csv(filename)

df = data_option.rename(columns={
    'indexprice': 'index_price',
    'strike_price': 'strike',
    'iv': 'mark_iv',
    'option_type': 'is_call',
    })
df['date'] = pd.to_datetime(df['date'])
df['is_call'] = np.where(df['is_call'] == 'C', 1, 0)

# Unique dates and tau values
unique_dates = df['date'].unique()
unique_taus = df['tau'].unique()  # Assuming 'tau' column exists in df

int_rate = 0

# keep_date = ["2022-11-17","2022-10-29","2022-07-22","2022-08-31"]
# filtered_df = df[df['date'].isin(pd.to_datetime(keep_date))]
# filtered_df.to_csv("CryptoRaw_data_2.csv", index=False)


# Create output directory
output_dir = 'Q_density_files_new'
os.makedirs(output_dir, exist_ok=True)


#######################################################
#
# Step 1: generating the Q density (original)
#
#######################################################

# Loop through each date and tau
start_date = np.datetime64('2021-03-29')
end_date = np.datetime64('2021-03-31')

unique_dates = unique_dates[(unique_dates>start_date) & (unique_dates<end_date)]
unique_dates = np.sort(unique_dates)
unique_taus = np.sort(unique_taus)

unique_dates = ['2022-11-17', '2022-10-29']
unique_dates = np.array(unique_dates, dtype='datetime64')

# date = np.datetime64('2022-09-08')

for date in unique_dates:
    sub_pre = df[(df['is_call'] == 1) & 
              (df['date'] == date) & 
               (df['moneyness'] >= 0.8) & (df['moneyness'] < 1.2) &
               (df['mark_iv'] > 0)]
       
    unique_taus_sub = sub_pre['tau'].unique()
    unique_taus_sub = np.sort(unique_taus_sub)
    date_str = np.datetime_as_string(date)[:10]
    
    for tau in unique_taus_sub:
        print(f"Processing date {date_str} and tau {tau}...")
        sub = sub_pre[(sub_pre['tau'] == tau)]
        if max(sub['moneyness']) > 1.10 and min(sub['moneyness']) < 0.85:
      
            try:
                # spd, sub = spdbl(df, date, tau, int_rate)
                spd, sub,spd_band = spdbl(sub, date, tau, int_rate)
            except ValueError as e:
                print(f"Skipping date {date_str} and tau {tau} due to error: {e}")
                continue
            
            # Format date and tau for filename
            # date_str = date.strftime('%Y%m%d')
           
            # tau_str = str(tau).zfill(3)  # Zero-fill tau to 3 digits if necessary
            tau_str = str(tau)
            if isinstance(spd, int):
                print('no result')
            else:
                output_filename = os.path.join(output_dir, f'raw_Q_density_{date_str}_tau{tau_str}.csv')
                spd.to_csv(output_filename, index=False)
                
            if isinstance(spd_band, int):
                print('no result')
            else:
                output_filename = os.path.join(output_dir, f'raw_Q_density_{date_str}_tau{tau_str}_band.csv')
                spd_band.to_csv(output_filename, index=False)




# # one example
# date = pd.to_datetime('2018-09-01')
# tau = 27
# int_rate = 0
# spd, sub  = spdbl(df, date, tau, int_rate)
# output_filename = 'Q_density_20180901.csv'
# spd.to_csv(output_filename, index=False)

# print(spd)
# spd.head()
# spd['y'].describe()

#######################################################
#
# Step 2: generating the Q density (tail fitting)
#
#######################################################

import pandas as pd
import numpy as np
from scipy.interpolate import splev, splrep
from scipy.stats import genextreme as gev

# path = '/Users/tracy/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/PAPER/Pricing/code/test/'
path = "C:/Users/leizhou/OneDrive - National University of Singapore/PAPER/Pricing/code/test/"
os.chdir(path)



# Define the Q_tail_logret_Figlewski function
def Q_tail_logret_Figlewski(spdy, m):
    # data filtering
    spdy = spdy[(m > 0.5) & (m <= 2)]
    m = m[(m > 0.5) & (m <= 2)]

    # Check monotonicity and filter data, Find the maximum spdy
    i = np.argmax(spdy)
    spdy1 = spdy[:i+1]
    spdy2 = spdy[i:]
    m1 = m[:i+1]
    m2 = m[i:]
    a1 = np.diff(spdy1)
    a2 = np.diff(spdy2)

    # get the final spdy and m
    spdy = np.concatenate((spdy1[np.concatenate(([True], a1 > 0))], spdy2[np.concatenate((a2 < 0, [True]))]))
    m = np.concatenate((m1[np.concatenate(([True], a1 > 0))], m2[np.concatenate((a2 < 0, [True]))]))
    
    r_t = np.log(m)
    q_rt = spdy * m

    def fit_tail(target, r_t, q_rt):
        k = min(3, len(r_t) - 1)  # Ensure k is less than the number of data points
        rnd = splev(target, splrep(r_t, q_rt, k=k))
        bounds = [(-0.1, 0.5), (0.01, 0.25), (-0.5, 0.5)]

        def fitness_function(x):
            return sum([
                abs(gev.pdf(target[0], x[0], loc=x[2], scale=x[1]) - rnd[0]),
                abs(gev.pdf(target[1], x[0], loc=x[2], scale=x[1]) - rnd[1]),
                abs(gev.pdf(target[2], x[0], loc=x[2], scale=x[1]) - rnd[2])
            ])

        result = differential_evolution(fitness_function, bounds)
        return result.x

    if np.max(m) <= 1.2:
        target_r = [0.18, 0.19, 0.2]
    elif np.max(m) <= 2:
        target_r = [np.log(np.max(m)) - 0.02, np.log(np.max(m)) - 0.01, np.log(np.max(m))]
    else:
        target_r = [0.98, 0.99, 1]

    solution_r = fit_tail(target_r, r_t, q_rt)

    if np.min(m) >= 0.8:
        target_l = [-0.22, -0.21, -0.2]
    elif np.min(m) >= 0:
        target_l = [np.log(np.min(m)), np.log(np.min(m)) + 0.01, np.log(np.min(m)) + 0.02]
    else:
        raise ValueError("Minimal moneyness is negative.")

    solution_l = fit_tail(target_l, r_t, q_rt)

    return_range = np.linspace(-1, 1, 200)
    q_r = gev.pdf(np.linspace(r_t[-1] + 0.001, max(return_range), 100), solution_r[0], loc=solution_r[2], scale=solution_r[1])
    q_l = gev.pdf(np.linspace(-min(return_range), -r_t[0] + 0.001, 100), solution_l[0], loc=solution_l[2], scale=solution_l[1])

    rt = np.concatenate((np.linspace(min(return_range), r_t[0] - 0.001, 100), r_t, np.linspace(r_t[-1] + 0.001, max(return_range), 100)))
    Q_rt = np.concatenate((q_l, q_rt, q_r))

    return rt, Q_rt


# List all files in the folder
folder_path = folder_path = os.path.join(path, "Q_density_files/")
files = os.listdir(folder_path)
output_files = [file for file in files if file.startswith('raw_Q')]


# Loop through each file in the folder
for file in output_files:
    input_filename = os.path.join(folder_path, file)
    data = pd.read_csv(input_filename)

    # Extract spdy and m columns
    spdy = data['y'].values
    m = data['m'].values

    # Call the function to process data
    rt, Q_rt = Q_tail_logret_Figlewski(spdy, m)

    # Write Q density to CSV
    output = pd.DataFrame({'Return': rt, 'Q_density': Q_rt})
    
    # Create output filename by removing 'raw_' from the input filename
    output_filename = file.replace('raw_', '')
    output_path = os.path.join(folder_path, output_filename)
    output.to_csv(output_path, index=False)
    print(f"Processed and saved: {output_filename}")



from tqdm import tqdm
import numpy as np
from scipy.optimize import differential_evolution

# Loop through each file in the folder with progress bar
for i, file in enumerate(tqdm(output_files, desc="Processing files", unit="file")):
    input_filename = os.path.join(folder_path, file)
    data = pd.read_csv(input_filename)

    # Extract spdy and m columns
    spdy = data['y'].values
    m = data['m'].values

    # Call the function to process data
    rt, Q_rt = Q_tail_logret_Figlewski(spdy, m)

    # Write Q density to CSV
    output = pd.DataFrame({'Return': rt, 'Q_density': Q_rt})
    
    # Create output filename by removing 'raw_' from the input filename
    output_filename = file.replace('raw_', '')
    output_path = os.path.join(folder_path, output_filename)
    output.to_csv(output_path, index=False)

    print(f"Processed and saved: {output_filename}")



