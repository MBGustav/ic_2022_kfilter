import numpy as np

# Parâmetros do sistema
F = 1.2
B = 0.5
H = 1.0

# Ruído do processo e de medição
process_noise_std = 0.1
measurement_noise_std = 0.5

# Função de transição do sistema
def system_dynamics(x, u):
    return F * x + B * u + np.random.normal(0, process_noise_std)

# Geração de dados de treinamento
num_samples = 1000
states = np.zeros(num_samples)
controls = np.random.uniform(-1, 1, num_samples)
measurements = np.zeros(num_samples)

for k in range(1, num_samples):
    states[k] = system_dynamics(states[k-1], controls[k])
    measurements[k] = H * states[k] + np.random.normal(0, measurement_noise_std)

# Implementação da Rede Neural
input_size = 1
hidden_size = 10
output_size = 2
learning_rate = 0.01

# Inicialização dos pesos da rede
weights_input_hidden = np.random.randn(input_size, hidden_size)
weights_hidden_output = np.random.randn(hidden_size, output_size)

# Treinamento da Rede Neural
for epoch in range(50):
    for k in range(1, num_samples):
        # Forward Pass
        hidden_layer = np.maximum(0, np.dot(states[k-1], weights_input_hidden))
        neural_output = np.dot(hidden_layer, weights_hidden_output)

        # Cálculo do erro
        error = states[k] - neural_output

        # Backward Pass
        # Atualização dos pesos
        weights_hidden_output += learning_rate * hidden_layer.reshape(-1, 1) * error
        weights_input_hidden += learning_rate * np.outer(states[k-1], np.dot(weights_hidden_output, error))

# Filtro de Kalman
def kalman_filter(x_pred, P_pred, z):
    # Predição
    x_hat = F * x_pred
    P_hat = F**2 * P_pred + process_noise_std**2
    
    # Atualização
    K = P_hat * H / (H**2 * P_hat + measurement_noise_std**2)
    x_updated = x_hat + K * (z - H * x_hat)
    P_updated = (1 - K * H) * P_hat
    
    return x_updated, P_updated

# Aplicação do Filtro de Kalman usando saídas da Rede Neural
filtered_states = np.zeros(num_samples)
P = 1.0  # Estimativa inicial da covariância do estado

for k in range(1, num_samples):
    # Usando a saída da Rede Neural como entrada para o Filtro de Kalman
    hidden_layer = np.maximum(0, np.dot(states[k-1], weights_input_hidden))
    neural_output = np.dot(hidden_layer, weights_hidden_output)

    # Aplicação do Filtro de Kalman
    filtered_states[k], P = kalman_filter(neural_output, P, measurements[k])

# Plotagem dos resultados
import matplotlib.pyplot as plt

plt.plot(states, label='True States', linestyle='dashed')
plt.plot(measurements, label='Measurements', marker='o', linestyle='None', alpha=0.7)
plt.plot(filtered_states, label='Filtered States', marker='o', linestyle='None', alpha=0.7)
plt.legend()
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('Kalman Filter with Neural Network Prediction')
plt.show()
