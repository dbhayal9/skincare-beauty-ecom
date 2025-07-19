module.exports = function(webpackEnv) {
  return {
    devServer: {
      allowedHosts: "all",
      host: 'localhost',
      port: 3000
    }
  }
}
