# TCC — mRMR e inferência de conectividade neuronal

**Murilo Cassiavilani | Estatística — UFSCar**

Repositório de organização dos códigos utilizados durante o Trabalho de Graduação A (TGA) e B (TGB). O objetivo é facilitar a localização dos experimentos, acompanhar suas versões e permitir a reprodução dos resultados quando os dados e dependências estiverem disponíveis.

## Onde encontrar

| Etapa | Pasta | Conteúdo |
|---|---|---|
| TGA | [Fundamentos](TGA/01_fundamentos/) | Cadeias de Markov e estimadores |
| TGA | [Classificação](TGA/02_classificacao/) | Wisconsin e Colon Cancer |
| TGA | [Redes neuronais](TGA/03_redes_neuronais/) | GL, seleção mRMR e circuito reduzido |
| TGA | [Dados anatômicos](TGA/04_dados_anatomicos/) | Extração de matrizes W |
| TGB | [Regra de parada](TGB/01_regra_parada/) | Entregável B |
| TGB | [Grafos simulados](TGB/02_grafos_simulados/) | Cenários 1–3, pós-processamento e mapas |
| TGB | [MaleCNS](TGB/03_malecns/) | Transformações de pesos e redes N12/N25/N50 |
| TGB | [ROI](TGB/04_roi/) | Reconstrução e comparação de T |
| TGB | [Monte Carlo](TGB/05_monte_carlo/) | Figuras do capítulo 6 |

## Reprodução dos mapas dos cenários

1. Os checkpoints da bateria v2 devem estar disponíveis em `resultados_4_1_v2_COMPLETA/cenario_N/T_N/rep_N.rds`.
2. Execute [02_resumir_checkpoints.R](TGB/02_grafos_simulados/02_resumir_checkpoints.R) se desejar tabelas de métricas e identificação correta.
3. O script [03_heatmaps_frequencia.R](TGB/02_grafos_simulados/03_heatmaps_frequencia.R) lê diretamente os checkpoints, sem simular novamente, e gera PNGs dos três cenários. Configure a variável de ambiente `TCC_RESULTADOS_GRAFOS` com o caminho da pasta se necessário.

## Notas importantes

- Os scripts originais são preservados como registros históricos, com exceção da remoção de uma credencial literal de neuPrint. Use `NEUPRINT_TOKEN` no ambiente.
- A anatomia do MaleCNS é real; os disparos são simulados. As versões do GL e os critérios de seleção variam entre experimentos.
- Alguns scripts têm caminhos absolutos da VM ou do computador original. Ajuste-os antes de executar.
- Não publique tokens, arquivos locais de configuração, checkpoints volumosos ou dados sem licença adequada.
- A execução completa das baterias pode ser demorada. Confirme as configurações e utilize os pilotos quando disponíveis.
- A organização não implica que todos os scripts tenham sido executados ou validados neste repositório. Os arquivos de dados e checkpoints da VM não foram enviados.
- Consulte o [inventário](docs/INVENTARIO.md) para detalhes sobre os scripts e dependências pendentes.
